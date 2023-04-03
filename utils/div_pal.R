
div_pal <- function(x, by, lower, upper) {
  x.breaks <- c(min(x),
                seq(from = lower, to = upper, by = by), max(x))
  # Palette and colors
  col_div <- pal_div(n = length(x.breaks) - 1)[cut(x = x,
                                                   breaks = x.breaks, right = FALSE,
                                                   include.lowest = TRUE, labels = FALSE)]
  return(list("col" = col_div, "breaks" = x.breaks))
}
pal_div <- c(
  brewer.pal(n = 9, name = "Blues") %>% rev,
  brewer.pal(n = 9, name = "Oranges")) %>%
  colorRampPalette()