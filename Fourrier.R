#Gaussian Blur by Fourier Transformation
library("jpeg")
library("MASS")

#Image files to be used
img1 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Minami.jpg", native=TRUE)
img2 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Veneza.jpg", native=TRUE)
img3 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Outono.jpg", native=TRUE)
img4 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Rosto.jpg", native=TRUE)
img5 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Zebra.jpg", native=TRUE)
img6 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Sol.jpg", native=TRUE)
img7 <- readJPEG("C:/Users/Renato Rissi/Documents/PIV/Pisa.jpg", native=TRUE)



#Creating the gaussian matrix
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
    A <- img[1:nrow(img)-1,1:nrow(img)-1]
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

imgf <- conv(img7,0.1)
image(Re(gmask(img7,2)), col = gray((1:255)/255))