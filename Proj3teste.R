library("jpeg")
library("squash")

img1 <- "C:/Users/Renato Rissi/Documents/PIV/Minami.jpg"
img2 <- "C:/Users/Renato Rissi/Documents/PIV/Veneza.jpg"
img3 <- "C:/Users/Renato Rissi/Documents/PIV/Outono.jpg"
img4 <- "C:/Users/Renato Rissi/Documents/PIV/Rosto.jpg"
img5 <- "C:/Users/Renato Rissi/Documents/PIV/Zebra.jpg"
img6 <- "C:/Users/Renato Rissi/Documents/PIV/Sol.jpg"
img7 <- "C:/Users/Renato Rissi/Documents/PIV/Pisa.jpg"


sobel <- function(imagem){
    A <- readJPEG(imagem, native=TRUE)
    gx <- matrix(0,254,254)
    gy <- matrix(0,254,254)
    mg <- matrix(0,254,254)
    for(j in 2:255){
        for(i in 2:255){
            gx[i-1,j-1] <- (A[i-1,j+1]+2*A[i,j+1]+A[i+1,j+1])-(A[i-1,j-1]+2*A[i,j-1]+A[i+1,j-1])
            gy[i-1,j-1] <- (A[i-1,j-1]+2*A[i-1,j]+A[i-1,j+1])-(A[i+1,j-1]+2*A[i+1,j]+A[i+1,j+1])
            mg[i-1,j-1] <- sqrt(gx[i-1,j-1]**2+gy[i-1,j-1]**2)
        }
    }
    mg <- t(mg)
    mg <- mg/max(mg)
    return(mg)
}

line <- function(mm, c, a){
    mtx <- matrix(0,17,17)
    for(i in 1:17){
        j=round(i*mm+c)
        k=round(i*mm+c+1)
        l=round(i*mm+c-1)
        if(18-j<18){
            mtx[18-j,i]=a
        }
        else{
            #do nothing
        }
    }
    return(mtx)
}

quat <- line(1,0,2)+line(1,1,-1)+line(1,-1,-1)
hori <- line(0,9,2)+line(0,10,-1)+line(0,8,-1)
trin <- line(sqrt(3)/3,4,2)+line(sqrt(3)/3,5,-1)+line(sqrt(3)/3,3,-1)
vert <- t(hori)
sess <- t(trin)

apl <- function(mask, A){
    tmask <- nrow(mask)                     #tamanho da mascara
    tA <- nrow(A)                           #tamanho da imagem
    meio <- floor(tmask/2)                  #1+meio=centro da matriz
    inicio <- ceiling(tmask/2)        #coordenada de inicio (centro da mascara)
    fim <- tA-meio                          #coordenada do fim
    imagem <- matrix(0,tA-(2*meio),tA-(2*meio))
    auxiliar <- 0
    for(j in inicio:fim){                   #percorre altura da imagem
        for(i in inicio:fim){               #percorre comprimento da imagem
            for(y in 0:(2*meio)){           #percorre a altura da mascara
                for(x in 0:(2*meio)){       #percorre o comprimento da mascara
                    auxiliar <- auxiliar + (mask[(y+1),(x+1)]*A[(j-meio+y),(i-meio+x)])
                }
            }
            imagem[j-meio,i-meio] <- sqrt(auxiliar**2)
            auxiliar <- 0
        }
    }
    imagem <- imagem/max(imagem)
    return(imagem)
}

#sob <- sobel(img1)
#plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
#rasterImage(sob,0,0,1,1)

img <- apl(quat, sob)
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(img,0,0,1,1)

plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(readJPEG("C:/Users/Renato Rissi/Documents/PIV/Minami.jpg", native=TRUE),0,0,1,1)
-------------------------------------------------------------------------

#teste sobel

> x <- matrix(0,3,3)
> x[1,1] <- 1
> x[1,2] <- 2
> x[1,3] <- 1
> x[3,3] <- -1
> x[3,2] <- -2
> x[3,1] <- -1
y <- (-1)*t(x)

imgteste <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Minami_mini.jpg", native=TRUE)

emx <- apl(x,imgteste)
emy <- apl(y,imgteste)
teste <- sqrt((emx**2)+(emy**2))
teste <- t(teste/max(teste))
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(teste,0,0,1,1)



