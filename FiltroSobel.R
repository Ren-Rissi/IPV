library(funr)
library(jpeg)
library(squash)

diretorio <- funr::detect_home()

img1 <- "Minami.jpg"
img2 <- "Veneza.jpg"
img3 <- "Outono.jpg"
img4 <- "Rosto.jpg"
img5 <- "Zebra.jpg"
img6 <- "Sol.jpg"
img7 <- "Pisa.jpg"




sobel <- function(imagem) {
    A <- readJPEG(imagem, native = TRUE)
    gx <- matrix(0, 254, 254)
    gy <- matrix(0, 254, 254)
    mg <- matrix(0, 254, 254)
#    fi <- matrix((,),254,254)

    for (j in 2:255){
        for (i in 2:255){
            gx[i - 1, j - 1] <-
            (A[i - 1, j + 1] + 2 * A [i, j + 1] + A[i + 1, j + 1]) -
            (A[i - 1, j - 1] + 2 * A[i, j - 1] + A[i + 1, j - 1])
            gy[i - 1, j - 1] <-
            (A[i - 1, j - 1] + 2 * A[i - 1, j] + A[i - 1, j + 1]) -
            (A[i + 1, j - 1] + 2 * A[i + 1, j] + A[i + 1, j + 1])
            mg[i - 1, j - 1] <- sqrt(gx[i - 1, j - 1]**2 + gy[i - 1, j - 1]**2)
#            fi[i-1,j-1] <- (gx[i-1,j-1]/mg[i-1,j-1],gy[i-1,j-1]/mg[i-1,j-1])
        }
    }
    mg <- t(mg)
    mg <- mg / max(mg)
#    plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
#    rasterImage(mg,0,0,1,1)
    return(mg)
}





gauss <- function(imagem){
    mask <- matrix(c(1, 1, 1, 1, 16, 1, 1, 1, 1), 3, 3)
    A <- readJPEG(imagem, native=TRUE)
#    A <- readJPEG(img1, native=TRUE)
    gx <- matrix(0, 254, 254)
    mg <- matrix(0, 254, 254)
    aux <- 0

    for (j in 2:255){
        for (i in 2:255){
            for (y in 1:3){
                for (x in 1:3){
                    aux <- aux + mask[x, y] * A[i - 2 + x, j - 2 + y]
                }
            }
            mg[i - 1, j - 1] <- sqrt(aux**2)
            aux <- 0
        }
    }

    mg <- t(mg)
    mg <- 1 - (mg / max(mg))
    plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
    rasterImage(mg, 0.5, 0.5, 1, 1)
    return(mg)
}



#------------------------ Projeto 3 --------------------------#

line <- function(mm, c, a) {
    mtx <- matrix(0, 17, 17)
    for (i in 1:17){
        j <- round(i * mm + c)
        k <- round(i * mm + c + 1)
        l <- round(i * mm + c - 1)
        if (18 - j < 18) {
            mtx[18 - j, i] <- a
        }
        else {
            #do nothing
        }
    }
    return(mtx)
}

quat <- line(1, 0, 2) + line(1, 1, -1) + line(1, -1, -1)
hori <- line(0, 9, 2) + line(0, 10, -1) + line(0, 8, -1)
trin <- line(sqrt(3) / 3, 4, 2) + line(sqrt(3) / 3, 5, -1) +
line(sqrt(3) / 3, 3, -1)
vert <- t(hori)
sess <- t(trin)


apl <- function(mask, A) {
    tmask <- nrow(mask) #tamanho da mascara
    tA <- nrow(A) #tamanho da imagem
    meio <- floor(tmask / 2) #1+meio=centro da matriz
    inicio <- ceiling(tmask / 2) #coordenada de inicio (centro da mascara)
    fim <- tA - meio #coordenada do fim
    imagem <- matrix(0, tA - (2 * meio), tA - (2 * meio))
    auxiliar <- 0
    for (j in inicio:fim){ #percorre altura da imagem
        for (i in inicio:fim){ #percorre comprimento da imagem
            for (y in 0:(2 * meio)){ #percorre a altura da mascara
                for (x in 0:(2 * meio)){ #percorre o comprimento da mascara
                    auxiliar <- auxiliar +
                    (mask[x + 1, y + 1] * A[i - meio + x, j - meio + y])
                }
            }
            imagem[i - meio, j - meio] <- sqrt(auxiliar**2)
            auxiliar <- 0
        }
    }
    imagem <- imagem / max(imagem)
    return(imagem)
}

sob <- sobel(img1)

img <- apl(hori, sob)
writeJPEG(img, target = "Minami_out.jpg", quality = 1, bg = "white")

    plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
    rasterImage(img, 0, 0, 1, 1)
