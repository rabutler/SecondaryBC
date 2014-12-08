qgaus=function(a,b,x){

w=1:5
xxy=1:5

w[1]=0.2955242247
w[2]=0.2692667193
w[3]=0.2190863625
w[4]=0.1494513491
w[5]=0.0666713443

xxy[1]=0.1488743389
xxy[2]=0.4333953941
xxy[3]=0.6794095682
xxy[4]=0.8650633666
xxy[5]=0.9739065285


xm=0.5*(b+a)
xr=0.5*(b-a)
ss=0
for (j in 1:5){
dx=xr*xxy[j]
xx1=xm+dx
xx2=xm-dx
ss=ss+w[j]*(integrand(xx1,x)+integrand(xx2,x))    #integrand is the function to integrate
}
ss=xr*ss

ss
}
