#' @title A function to mix two different colors
#'
#' @description This function mixes two different colors. It relies on functions of the package colorspace
#' @param alpha value between 0 and 1 where 0 is 100% the first color and 1 is 100% the second color.
#' @param color1 First color to be mixed
#' @param color2 Second color to be mixed

#'
#' @export
#'

mix_color <- function(alpha,color1,color2){
  new_color <- colorspace::mixcolor(alpha = alpha,
                                    color1=colorspace::RGB(col2rgb(color1)[1,1], col2rgb(color1)[2,1], col2rgb(color1)[3,1]),
                                    color2=colorspace::RGB(col2rgb(color2)[1,1], col2rgb(color2)[2,1], col2rgb(color2)[3,1]))
  new_color <- rgb(new_color@coords[1,1],new_color@coords[1,2],new_color@coords[1,3], alpha = 255, maxColorValue = 255)
  return(new_color)
}



