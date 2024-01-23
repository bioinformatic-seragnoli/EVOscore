format_segments <- function(seg){
  
  Fseg <- seg
  
  # format chr column
  chr_detect <- names(seg) %>% str_detect("[Cc]hr")
  if(sum(chr_detect)==1){
    idx <- names(seg) %>% str_detect("[Cc]hr") %>% which()
    oldname <- names(seg)[idx]
    Fseg <- Fseg %>% dplyr::rename("chr" = !!oldname)
  }
  
  # format start column
  start_detect <- names(seg) %>% str_detect("[Ss]tart")
  if(sum(start_detect)==1){
    idx <- names(seg) %>% str_detect("[Ss]tart") %>% which()
    oldname <- names(seg)[idx]
    Fseg <- Fseg %>% dplyr::rename("start" = !!oldname)
  }
  
  # format end column
  end_detect <- names(seg) %>% str_detect("[Ee]nd")
  if(sum(end_detect)==1){
    idx <- names(seg) %>% str_detect("[Ee]nd") %>% which()
    oldname <- names(seg)[idx]
    Fseg <- Fseg %>% dplyr::rename("end" = !!oldname)
  }
  
  # format n_probes column
  pro_detect <- names(seg) %>% str_detect("[Pp]robe|[Mm]arker")
  if(sum(pro_detect)==1){
    idx <- names(seg) %>% str_detect("[Pp]robe|[Mm]arker") %>% which()
    oldname <- names(seg)[idx]
    Fseg <- Fseg %>% dplyr::rename("n_probes" = !!oldname)
  }
  
  Fseg
  
}
