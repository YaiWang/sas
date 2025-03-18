library(dplyr)
library(lubridate)
library(purrr)
library(openxlsx)
library(here)
library(sf)
library(spdep)
library(ggplot2)
library(tmap)
library(rgdal)
library(ggplot2)
library(maptools)
library(mapproj)
library(CARBayes)
library(CARBayesST)
library(mgcv)  
library(geosphere)
library(proxy)

brazil_shape <- here("E:\\A.Research topic\\Dengue\\Maps\\bra_admbnda_adm1_ibge_2020.shp") %>%
  st_read()
# 给地图数据添加州名缩写 进行排序 方便之后和病例数据匹配
brazil_shape$ADM1_PT_1 <- c("AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", "MT", "MS",
                            "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO")
brazil_shape <- brazil_shape[order(brazil_shape$ADM1_PT_1), ]

filter_vec <- c("MG","RJ","SP","PR","RS","SC","GO","MT","MS","DF")

# 使用筛选条件进行行筛选
brazil_shape1 <- brazil_shape[brazil_shape$ADM1_PT_1 %in% filter_vec, ]

PR_map <- here("E:\\A.Research topic\\Dengue\\research1\\maps\\PR-10\\PR_Mesorregioes_2022.dbf") %>%
  st_read()
# 提取质心坐标
centroids <- st_centroid(PR_map)
plot(PR_map[-c(8,9),])
plot(PR_map[c(2,3,6,7,8,9),])
nb_q <- poly2nb(st_geometry(PR_map), queen = TRUE)
# 提取质心坐标的经度和纬度
centroid_coords <- st_coordinates(centroids)

# 创建包含经度和纬度的新数据框
PR_l <- data.frame(
  region = PR_map$NM_MESO,
  Longitude = centroid_coords[, "X"],  # 提取经度
  Latitude = centroid_coords[, "Y"]    # 提取纬度
)


###### 计算距离 ######
# 经纬度
lon_lat <- data.frame(
  region = paste0("R",1:10),
  Longitude = brazil_shape1$COORD_X,
  Latitude = brazil_shape1$COORD_Y
)

lon_lat_zong <- rbind(PR_l,lon_lat)
lon_lat_zong <- lon_lat_zong[-16,]

euclidean_distance <- function(lon1, lat1, lon2, lat2) {
  return(sqrt((lon2 - lon1)^2 + (lat2 - lat1)^2))
}

coords <- lon_lat_zong[, c("Longitude", "Latitude")]
distance_matrix <- dist(coords)
distance_matrix1 <- distm(coords, fun = distHaversine)
# Distance in kilometers 球面距离
great_circle_distance <- function(lon1, lat1, lon2, lat2) {
  
  # Convert degrees to radians
  lon1 <- lon1 * (pi / 180)
  lat1 <- lat1 * (pi / 180)
  lon2 <- lon2 * (pi / 180)
  lat2 <- lat2 * (pi / 180)
  
  # Earth radius in kilometers
  R <- 6371
  
  # Calculate differences
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  # Calculate great-circle distance using Haversine formula
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * asin(sqrt(a))
  
  # Distance in kilometers
  distance <- R * c
  
  return(distance)
}

# 创建一个空的距离矩阵
num_locations <- nrow(lon_lat_zong)
distance_matrix <- matrix(0, nrow = num_locations, ncol = num_locations)

# 计算每对地区之间的距离
for (i in 1:(num_locations - 1)) {
  for (j in (i + 1):num_locations) {
    # 获取地区i和地区j的经纬度
    lon1 <- lon_lat_zong$Longitude[i]
    lat1 <- lon_lat_zong$Latitude[i]
    lon2 <- lon_lat_zong$Longitude[j]
    lat2 <- lon_lat_zong$Latitude[j]
    
    # 计算距离
    distance <- great_circle_distance(lon1, lat1, lon2, lat2)
    
    # 将距离存储在距离矩阵中
    distance_matrix[i, j] <- distance
    distance_matrix[j, i] <- distance
  }
}


# 人口数据
pop <- readxl::read_excel("E:/A.Research topic/Dengue/data/dengue_cases/Population_data.xlsx")
pop1 <- subset(pop, select = filter_vec)

pop2 <- pop1[order(names(pop1))]
N_h <- pop2[469:520,]
# PR pop 
pr_pop <- readxl::read_excel("E:\\A.Research topic\\Dengue\\research2\\code\\PR-10-POP.xlsx")
pop_s <- cbind(t(pr_pop$pop),N_h[1,])
colnames(pop_s) <- c(t(pr_pop$M1),colnames(N_h))
pop_s <- pop_s[,-16]



# 流动总人口
T_i <- ceiling(pop_s * 0.1)
T_i <- data.matrix(T_i)
###### 辐射模型 ######

# 计算s_ij人口
pop_s <- data.matrix(pop_s)
population_in_circle_matrix <- matrix(0, nrow = 19, ncol = 19)

# 遍历每一对地区（i 和 j）
for (i in 1:19) {
  for (j in 1:19) {
    if (i != j) {
      distance <- distance_matrix[i, j]
      total_population <- 0
      
      # 遍历每个地区 k
      for (k in 1:19) {
        # 如果 k 不是 i 和 j，并且与 i 的质心距离小于或等于 distance，则将 k 的人口加入总和
        if (k != i && k != j && distance_matrix[i, k] <= distance) {
          total_population <- total_population + pop_s[k]
        }
      }
      
      # 将结果填充到矩阵中
      population_in_circle_matrix[i, j] <- total_population
    }
  }
}

s_ij <- data.matrix(population_in_circle_matrix)

T_ij <- matrix(0, nrow = 19, ncol = 19)

for (i in 1:19) {
  for (j in 1:19) {
    if (i != j){
      T_ij[i,j] <- ceiling(T_i[i]*(pop_s[i]*pop_s[j])/((pop_s[i]+s_ij[i,j])*(pop_s[i]+pop_s[j]+s_ij[i,j])))
    }
  }
}
setwd("E:\\A.Research topic\\Dengue\\research2\\code")
save(T_ij, file = "flow_pop.RData")




PR_cases <- readxl::read_excel("PR_cases.xlsx")

nR <- ncol(PR_cases)
nT <- nrow(PR_cases)
load("temp_quan.RData")

PR_temp <- temp$PR_TEMP
PR_temp1 <- matrix(PR_temp, nrow=nT, ncol=nR, byrow = FALSE)
temp_s <- cbind(PR_temp1,temp)
# This is the general function for the Briere fit.
briere <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*x*(x-T0)*sqrt(Tm-x)
}

# This is the general function for the quadratic fit. 
quadratic <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    0.0
  else
    c*(x-T0)*(x-Tm)
}

# This is the general function for the inverted quadratic fit.
inverted_quadratic <- function(x, c, T0, Tm){
  if((x < T0) | (x > Tm))
    1.0/24
  else
    1.0/(c*(x-T0)*(x-Tm))
}
# 叮咬率
a_f <- function(temp){
  briere(temp,2.02e-04,13.35,40.08)
}
# 蚊子感染概率
c_f <- function(temp){
  briere(temp,4.91e-04,12.22,37.46)
}
# 蚊子死亡率
mu_v_f <- function(temp){
  inverted_quadratic(temp,-1.48e-01, 9.16, 37.73)
}

a <- data.frame(matrix(NA, ncol = nR, nrow = nT))
c <- data.frame(matrix(NA, ncol = nR, nrow = nT))
mu_v <- data.frame(matrix(NA, ncol = nR, nrow = nT))

for(i in 1:nR){
  for(t in 1:nT){
    # biting rate
    a[t,i] <- 7*a_f(temp_s[t,i])
    
    # probability	of mosquito	infection per	bite	on	an	infectious	host
    c[t,i] <- c_f(temp_s[t,i])
    
    # the mosquito mortality rate
    mu_v[t,i] <- 7*mu_v_f(temp_s[t,i])
    
  }  
}

N_h <- matrix(pop_s, nrow=nT, ncol=nR, byrow = TRUE)

birth_rate <- readxl::read_excel("E:\\A.Research topic\\Dengue\\data\\dengue_cases\\birth rate.xlsx")

birth_rate1 <- birth_rate[rep(1:nrow(birth_rate),each=52),]
# 添加名为"week"的列，数据为1到52，重复十次
birth_rate1$week <- rep(1:52, times = 10)
birth_rate1 <- birth_rate1[, c(1, ncol(birth_rate1), 2:(ncol(birth_rate1)-1))]
birth_rate2 <- birth_rate1[469:520,]
rho <- birth_rate2$rho
mu_h <- birth_rate2$mu_h

rho <- rho/1000/365*7
mu_h <- mu_h/1000/365*7

save(T_ij, PR_cases, N_h, a, c, mu_v, rho, mu_h, nR, nT, file = "dengue_1.RData")
