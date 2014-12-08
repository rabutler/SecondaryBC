#"zfill" <- function(x,width){
"zfill" <- function(x){
  #x - integer input; width - width of return string with pre-padded zeros
  if (x < 10) {s=paste("00",x,sep="");return(s)}
  if ((x>=10)&&(x<100)) {s=paste("0",x,sep="");return(s)}
  return(x)


}