###############################
### Fisheries data analysis ###
###############################

# Landings for LFAs 34-38
landings <- read_csv("Landings34-38.csv", col_select = c(YR, LAND))
ggplot(filter(landings, YR %in% 2005:2023), aes(YR, LAND))+
  geom_col()+
  ylim(0, 50000)+
  labs(x = "Year", y = "Landings (kt)",
       title = "Lobster landings in LFAs 34-38")+
  theme_light()

# Now look at fisheries data within the subareas
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

landings <- readRDS("D:/Projects/Lobster/Landings_byGrid_year_34 1.rds")
grid <- readRDS("D:/Projects/Lobster/GridPolys_DepthPruned_37Split.rds")
grid <- grid[grid$LFA == 34,]
landings$GRID_NO <- as.character(landings$GRID_NO)

# Get the coordinates into the landings object
landings.sf <- left_join(landings,grid,by=c("GRID_NO","LFA"))
landings.sf <- st_as_sf(landings.sf)

landings.sf <- landings.sf |> group_by(SYEAR) |> 
  mutate(tot.land = sum(Landings,na.rm=T),
         tot.hauls = sum(TrapHauls,na.rm=T),
         tot.cpue = sum(Landings,na.rm=T)/sum(TrapHauls,na.rm=T)) #|> as.data.frame()

ggplot(landings.sf) + geom_sf(aes(fill = sqrt(Landings))) + facet_wrap(~SYEAR) + scale_fill_viridis_b()

ggplot(landings.sf) + geom_sf(aes(color = GRID_NO)) + geom_sf_text(aes(label = GRID_NO)) 

# These are the grid numbers that define the subareas
# SMB                 = 69,81,92
# Outer SMB           = 103, 114, 125, 126
# Greater Lobster Bay = 127, 140,141, 156, 157
# Jacquards Ridge     =  138,139,154,155

sub.areas <- as.character(c(69,81,92,103,114,125:127,138:141,154:157))
landings.sf$box <- "Outside"
landings.sf$box[landings.sf$GRID_NO %in% sub.areas] <- "Inside"

boxes.in.out <- landings.sf |> collapse::fgroup_by(box,SYEAR,tot.land,tot.hauls,tot.cpue) |>
  collapse::fsummarise(land.in.out = sum(Landings,na.rm=T),
                       hauls.in.out = sum(TrapHauls,na.rm=T),
                       cpue.in.out =  sum(Landings,na.rm=T)/sum(TrapHauls,na.rm=T)) |> as.data.frame()

boxes.in.out$prop.land <- boxes.in.out$land.in.out / boxes.in.out$tot.land
boxes.in.out$prop.hauls <- boxes.in.out$hauls.in.out / boxes.in.out$tot.hauls
boxes.in.out$per.cpue <- 100*((boxes.in.out$cpue.in.out - boxes.in.out$tot.cpue) /  boxes.in.out$tot.cpue)

ggplot(boxes.in.out) + geom_line(aes(x=SYEAR,y=land.in.out,color=box),lwd=2)

ggplot(boxes.in.out) + geom_line(aes(x=SYEAR,y=prop.land,color=box),lwd=2)
ggplot(boxes.in.out) + geom_line(aes(x=SYEAR,y=prop.hauls,color=box),lwd=2)
ggplot(boxes.in.out) + geom_line(aes(x=SYEAR,y=per.cpue,color=box),lwd=2)

# We make an object with the subareas in it
#landings.sub <- landings.sf |> subset(GRID_NO %in% sub.areas) #& SYEAR %in% 2005:2023)
landings.sf$name <- "34 Outside"
landings.sf$name[landings.sf$GRID_NO %in% c("69", "81", "92")] <- "SMB"
landings.sf$name[landings.sf$GRID_NO %in% c("103", "114", "125", "126")] <- "Outer SMB"
landings.sf$name[landings.sf$GRID_NO %in% c("127", "140","141", "156", "157")] <- "GLB"
landings.sf$name[landings.sf$GRID_NO %in% c("138","139","154","155")] <- "JR"

landings.sf$name <- factor(landings.sf$name,levels = c("SMB","Outer SMB","JR","GLB","34 Outside"))

landings.sub <- landings.sf |> collapse::fgroup_by(name,SYEAR,LFA,tot.land,tot.hauls,tot.cpue) |> 
  collapse::fsummarise(area.land = sum(Landings,na.rm=T),
                       area.hauls = sum(TrapHauls,na.rm=T),
                       area.cpue = sum(Landings,na.rm=T)/sum(TrapHauls,na.rm=T)) |> as.data.frame()

# Plots of landings
land.p <- ggplot(landings.sub[landings.sub$SYEAR <= 2023,]) + geom_line(aes(x=SYEAR,y=area.land/1e6,color=name),linewidth = 1.5) + 
  scale_color_manual(values = c("blue","black","grey","orange",'firebrick2')) + 
  scale_x_continuous(name="") + scale_y_continuous(name = "Landings (tonnes x 1000)",limits=c(0,NA)) +
  theme_bw(base_size = 18) + theme(legend.position = 'none') 



land.bp <- ggplot(landings.sub[landings.sub$SYEAR <= 2023,]) +  geom_bar(stat = 'identity',aes(x=SYEAR,y=area.land,fill=name),position='fill') + 
  scale_fill_manual(values = c("blue","black","grey","orange",'firebrick2')) + 
  scale_x_continuous(name="") + scale_y_continuous(name = "Proportion of Landings ") +
  theme_bw(base_size = 18) + theme(legend.title = element_blank()) 

land.p.combo <- plot_grid(land.p,land.bp,nrow=1,rel_widths = c(0.65,1))

save_plot("D:/Projects/Lobster/LFA 34 landings by subarea.png",land.p.combo,base_height = 8,base_width = 11)

# Plots of hauls
haul.p <- ggplot(landings.sub[landings.sub$SYEAR <= 2023,]) + geom_line(aes(x=SYEAR,y=area.hauls,color=name),linewidth = 1.5) + 
  scale_color_manual(values = c("blue","black","grey","orange",'firebrick2')) +
  scale_x_continuous(name="") + scale_y_continuous(name = "Effort (Number of hauls)",limits=c(0,NA)) +
  theme_bw(base_size = 18) + theme(legend.position = 'none') 

haul.bp <- ggplot(landings.sub[landings.sub$SYEAR <= 2023,]) +  geom_bar(stat = 'identity',aes(x=SYEAR,y=area.hauls,fill=name),position='fill') + 
  scale_fill_manual(values = c("blue","black","grey","orange",'firebrick2')) + 
  scale_x_continuous(name="") + scale_y_continuous(name = "Proportion of Effort ") +
  theme_bw(base_size = 18) + theme(legend.title = element_blank()) 


haul.p.combo <- plot_grid(haul.p,haul.bp,nrow=1,rel_widths = c(0.65,1))

save_plot("D:/Projects/Lobster/LFA 34 hauls by subarea.png",haul.p.combo,base_height = 8,base_width = 11)

# Finally the CPUE, it's just a single panel as proportion doesn't make sense here.
cpue.p <- ggplot(landings.sub[landings.sub$SYEAR <= 2023,]) + geom_line(aes(x=SYEAR,y=area.cpue,color=name),linewidth = 1.5) + 
  scale_color_manual(values = c("blue","black","grey","orange",'firebrick2'))+
  scale_x_continuous(name="") + scale_y_continuous(name = "CPUE",limits=c(0,NA)) +
  theme_bw(base_size = 18) + theme(legend.title = element_blank()) 

save_plot("D:/Projects/Lobster/LFA 34 cpue by subarea.png",cpue.p,base_height = 8,base_width = 8)

# What is the change from the start and end of the time series?
landings.three <- landings.sub |> collapse::fsubset(SYEAR %in% c(1999,2005,2023))


landings.three$prop.of.total.land <- landings.three$area.land/landings.three$tot.land
landings.three$prop.of.total.hauls <- landings.three$area.hauls/landings.three$tot.hauls
landings.three$per.of.average.cpue <- 100*((landings.three$area.cpue-landings.three$tot.cpue)/landings.three$tot.cpue)

# How much of the effort came out of these four areas
prop.land.of.LFA.2005 <- sum(landings.three$prop.of.total.land[landings.three$SYEAR == 2005 & landings.three$name != "34 Outside"])
prop.land.of.LFA.2023 <- sum(landings.three$prop.of.total.land[landings.three$SYEAR == 2023 & landings.three$name != "34 Outside"])
# Hauls
prop.effort.of.LFA.2005 <- sum(landings.three$prop.of.total.hauls[landings.three$SYEAR == 2005 & landings.three$name != "34 Outside"])
prop.effort.of.LFA.2023 <- sum(landings.three$prop.of.total.hauls[landings.three$SYEAR == 2023 & landings.three$name != "34 Outside"])


# Now to get the proportional changes
prop.change <- landings.three %>% select(name,SYEAR,area.land,area.hauls,area.cpue) %>%
  filter(SYEAR %in% c(2005, 2023)) %>%
  pivot_wider(names_from = SYEAR,
              values_from = c(area.land, area.hauls, area.cpue),
              values_fill = NA,names_sep =".")


prop.change$per.change.land  <-  100*((prop.change$area.land.2023 - prop.change$area.land.2005)/ prop.change$area.land.2005)
prop.change$per.change.hauls =  100*((prop.change$area.hauls.2023   - prop.change$area.hauls.2005)  / prop.change$area.hauls.2005)
prop.change$per.change.cpue  = 100*((prop.change$area.cpue.2023   - prop.change$area.cpue.2005)  / prop.change$area.cpue.2005)

# Proportional change in landings in the 34 Outside between 1999 and the peak year of 2016
100*((landings.sub$area.land[landings.sub$SYEAR == 2016 & landings.sub$name == "34 Outside"]/
        landings.sub$area.land[landings.sub$SYEAR == 1999 & landings.sub$name == "34 Outside"]) -1)
# The decline from 2016 to 2023 was...
100*((landings.sub$area.land[landings.sub$SYEAR == 2023 & landings.sub$name == "34 Outside"]/
        landings.sub$area.land[landings.sub$SYEAR == 2016 & landings.sub$name == "34 Outside"]) -1)
# Now each of the other 4 areas, these all declined over time, so just getting change from 1999 to 2023
100*((landings.sub$area.land[landings.sub$SYEAR == 2023 & landings.sub$name == "SMB"]/
        landings.sub$area.land[landings.sub$SYEAR == 1999 & landings.sub$name == "SMB"]) -1)
100*((landings.sub$area.land[landings.sub$SYEAR == 2023 & landings.sub$name == "Outer SMB"]/
        landings.sub$area.land[landings.sub$SYEAR == 1999 & landings.sub$name == "Outer SMB"]) -1)
100*((landings.sub$area.land[landings.sub$SYEAR == 2023 & landings.sub$name == "JR"]/
        landings.sub$area.land[landings.sub$SYEAR == 1999 & landings.sub$name == "JR"]) -1)
100*((landings.sub$area.land[landings.sub$SYEAR == 2023 & landings.sub$name == "GLB"]/
        landings.sub$area.land[landings.sub$SYEAR == 1999 & landings.sub$name == "GLB"]) -1)

# What did CPUE do, this is from 1999 to peak CPUE year, so year will differ by area.
# The decline from 2016 to 2023 was...
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2016 & landings.sub$name == "34 Outside"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 1999 & landings.sub$name == "34 Outside"]) -1)
# Now each of the other 4 areas, these all declined over time, so just getting change from 1999 to 2023
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2017 & landings.sub$name == "SMB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 1999 & landings.sub$name == "SMB"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2014 & landings.sub$name == "Outer SMB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 1999 & landings.sub$name == "Outer SMB"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2015 & landings.sub$name == "JR"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 1999 & landings.sub$name == "JR"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2016 & landings.sub$name == "GLB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 1999 & landings.sub$name == "GLB"]) -1)
# Subsequent declines...
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2023 & landings.sub$name == "34 Outside"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 2016 & landings.sub$name == "34 Outside"]) -1)
# Now each of the other 4 areas, these all declined over time, so just getting change from 1999 to 2023
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2023 & landings.sub$name == "SMB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 2017 & landings.sub$name == "SMB"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2023 & landings.sub$name == "Outer SMB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 2014 & landings.sub$name == "Outer SMB"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2023 & landings.sub$name == "JR"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 2015 & landings.sub$name == "JR"]) -1)
100*((landings.sub$area.cpue[landings.sub$SYEAR == 2023 & landings.sub$name == "GLB"]/
        landings.sub$area.cpue[landings.sub$SYEAR == 2016 & landings.sub$name == "GLB"]) -1)

View(prop.change)