library("jpeg")
library("squash")

#---------------------------LIMIARIZAÇÃO----------------------------------

#Recebe a imgem
#Verifica as dimensões
#Cria uma matriz nova de mesmo tamanho
#Separa o fundo do objeto e põe na nova matriz

lm <- function(img){
   tam <- dim(img)
   li <- matrix(0, tam[1], tam[2])
   for(i in 1:tam[1]){
      for(j in 1:tam[2]){
         if(img[i,j] <= 0.68) li[i,j] = 1
      }
   }
#image(lim, col = gray((0:255)/255))
return(li)
}

#--------------------------------PCA-------------------------------------

#alinha as folhas
#problema com prcomp
#modo bruteforce?
#(matriz covariancia, autovalor, autovetor, ordenação, reconstrução)

#-------------------------------Área-------------------------------------

#Recebe imagem limiarizada
#Conta os pixels que fazem parte do objeto como medida de área
#Se o pixel faz parte do objeto, contador +1

ar <- function(lim){
   area <- 0
   tam <- dim(lim)
   for(i in 1:tam[1]){
      for(j in 1:tam[2]){
      if(lim[i,j] > 0){area = area+1}
      }
   }
   return(area)
}

#------------------------------Perímetro 4----------------------------------
#Recebe a imagem limiarizada
#Testa os pixels da imagem pra descobrir se é borda
#Cria uma lista de pixels que são borda
#Conta os itens da lista

pr4 <- function(lim){
   tam <- dim(lim)
   per4 <- 0
   for(i in 2:(tam[1]-1)){
      for(j in 2:(tam[2]-1)){
         if(lim[i,j]>0){ if(lim[i-1,j]==0 ||
                            lim[i+1,j]==0 ||
                            lim[i,j-1]==0 ||
                            lim[i,j+1]==0){per4 <-  per4+1 } }
      }
   }
   return(per4)
}

#-----------------------------Perímetro 8-----------------------------------

pr8 <- function(lim){
   tam <- dim(lim)
   per8 <- 0
   for(i in 2:(tam[1]-1)){
      for(j in 2:(tam[2]-1)){
         if(lim[i,j]>0){ if(lim[i-1,j-1]==0 ||
                            lim[i-1,j]==0 ||
                            lim[i-1,j+1]==0 ||
                            lim[i,j-1]==0 ||
                            lim[i,j+1]==0 ||
                            lim[i+1,j-1]==0 ||
                            lim[i+1,j]==0 ||
                            lim[i+1,j+1]==0){per8 <-  per8 +1 } }
      }
   }
   return(per8)
}

#-----------------------------Diâmetro--------------------------------------
#Recebe a imagem limiarizada e rotacionada por PCA
#Encontra o primeiro e último pixel do objeto
#calcula a distância de tais pixels

dA <- function(lim){
   tam <- dim(lim)
   xi <- 0
   yi <- 0

   for(i in 1:tam[1]){
      for(j in 1:tam[2]){
         if(lim[i,j]>0){
            xi=i
            yi=j
            break
         }
      }
   }

   xf <- 0
   yf <- 0

   for(i in tam[1]:1){
      for(j in tam[2]:1){
         if(lim[i,j]>0){
            xf <- i
            yf <- j
            break
         }
      }
   }

   diamA <- sqrt(((xf-xi)**2)+((yf-yi)**2))
   return(diamA)
}

#-------------------------------------------------Outra direção

dB <-function(lim){
   tam <- dim(lim)
   xi <- 0
   yi <- 0

   for(i in 1:tam[2]){
      for(j in 1:tam[1]){
         if(lim[j,i]>0){
            xi=i
            yi=j
            break
         }
      }
   }

   xf <- 0
   yf <- 0

   for(i in tam[2]:1){
      for(j in tam[1]:1){
         if(lim[j,i]>0){
            xf <- i
            yf <- j
         break
         }
      }
   }

   diamB <- sqrt(((xf-xi)**2)+((yf-yi)**2))
   return(diamB)
}
#-----------------------------------------ecentricidade---------------------
#tomando de o PCA vai rotacionar a imagem de modo que a maior variação fique
#no eixo X, usamos a medida dos diâmetros para calcular a excentricidade

ec <- function(da,db){
   da -> diamA
   db -> diamB
   e<- sqrt(1-((diamB**2)/(diamA**2)))
   return(e)
}

#----------------------------Entropia---------------------------------------
entropy <- function(img){
   hit <- matrix(0, 256)
   x <- matrix(0:255)
   d <- dim(img)
   t <- length(img)
   for(i in 1:d[1]){
      for(j in 1:d[2]){
         hit[floor(img[i,j]*255)+1] <- hit[floor(img[i,j]*255)+1]+1
      }
   }

   p <- hit/t
   ent <- 0
   aux <- 0
   for(i in 1:256){
      if(p[i]!=0){ent <- ent + p[i]*log2(1/p[i])}
   }
   return(ent)
}
#---------------------------Filtro gaussiana borrar ------------------------
gmask <- function(image,sigma){
    nco <- ncol(image)-1
    midd <- ceiling(nco/2)
    mask <- matrix(1, nco, nco)
    for(y in 1:nco){
        for (x in 1:nco){
            xx=(x-(midd))+0.5
            yy=(y-(midd))+0.5
            mask[y,x]= exp(((xx**2+yy**2)*(-1))/(2*(sigma**2)))/(2*pi*(sigma**2))
        }
    }
    chess <- matrix(0,ncol=1,nrow=nco)
    for(z in 1:nco){chess[z]<-(-1)**(z+1)}
    chess <- (chess %*% t(chess))
    fmask <- fft(mask)*chess
    fmask <- Re(fft(fmask,inverse=TRUE)/length(fmask))
    return(abs(fmask))
}

conv <- function(img,sigma){
    A <- img[1:nrow(img)-1,1:nrow(img)-1]#####################################
    B <- gmask(img,sigma)
    fA <- fft(A)
    fB <- fft(B)
    C <- fA*fB
    C <- Re(fft(C, inverse=TRUE)/length(C))
    C <- abs(C)
    C <- t(C/max(C))
    plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
    rasterImage(C,0,0,1,1)
    return(C)
}


#---------------------------Tabela de medidas-------------------------------

tabela <- matrix(0, 60, 8)

for(i in 1:60){
   a <- "C:/Users/Renato Rissi/Documents/PIV/Folhas/img"
   b <- format(i)
   c <- ".jpg"
   img <- paste(a,b,c, sep = "", collapse = "")
   img <- readJPEG(img)
   lim <- lm(img)
   
   tabela[i,1] <- ar(lim)
   tabela[i,2] <- pr4(lim)
   tabela[i,3] <- pr8(lim)
   tabela[i,4] <- dA(lim)
   tabela[i,5] <- dB(lim)
   tabela[i,6] <- entropy(img)
#   tabela[i,7] <- entropy(conv(img,0.1))
#   tabela[i,8] <- entropy(conv(img,1))
}
image(img, col = gray((0:255)/255))

