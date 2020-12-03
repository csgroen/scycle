from plotnine import theme_light, theme, element_blank, element_text, element_rect

theme_std = theme_light() + theme(panel_grid=element_blank(), 
                                  strip_text = element_text(color = 'black'),
                                  strip_background=element_rect(fill = 'white'))
