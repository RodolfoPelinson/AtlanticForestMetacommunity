#' @title A function lightens a given color
#'
#' @description This function lightens a given color. It basically mixes a given color with white.
#' @param alpha value between 0 and 1 where 0 is white and 1 is 100% the chosen color.
#' @param color The color to be lightened

#'
#' @export
#'

lighten_color <- function(alpha,color){
  new_color <- colorspace::mixcolor(alpha = alpha,
                                    color1=colorspace::RGB(col2rgb(color)[1,1], col2rgb(color)[2,1], col2rgb(color)[3,1]),
                                    color2=colorspace::RGB(col2rgb("white")[1,1], col2rgb("white")[2,1], col2rgb("white")[3,1]))
  new_color <- rgb(new_color@coords[1,1],new_color@coords[1,2],new_color@coords[1,3], alpha = 255, maxColorValue = 255)
  return(new_color)
}
