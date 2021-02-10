# moorelr/Melts-Vols is licensed under The MIT License
# Copyright 2019 Lowell R. Moore


# Start ####

make_beep <- FALSE       # Make a beep noise at the end. Requires "beepr" package.
status_updates <- TRUE   # Print status updates while running.
pause_time <- 2          # Delay time if status updates are used

if(make_beep){library(beepr)}
print_status <- function(message){
  if(status_updates){
    print(message)
    Sys.sleep(pause_time)
  }
}

print_status("Setting things up...")

settings <- read.delim("SETTINGS.txt", sep = "=", header = FALSE, stringsAsFactors = FALSE)

# Default settings
norm_volume <- as.logical(settings[1, 2])
exclude_liquid <- as.logical(settings[2, 2])

get_info <- function(string){
  
  # get name
  start_pos <- 1
  stop_pos <- gregexpr(text = string, pattern = "   ")[[1]][1] - 1
  name <- substr(x = string, start = start_pos, stop = stop_pos)
  
  # get mass
  start_pos <- gregexpr(text = string, pattern = "mass = ")[[1]][1] + 7
  stop_pos <- gregexpr(text = string, pattern = "\\(gm\\)")[[1]][1] - 2
  mass <- substr(x = string, start = start_pos, stop = stop_pos)
  
  # get density
  start_pos <- gregexpr(text = string, pattern = "density = ")[[1]][1] + 10
  stop_pos <- gregexpr(text = string, pattern = "\\(gm/cc\\)")[[1]][1] - 2
  density <- substr(x = string, start = start_pos, stop = stop_pos)

  return(c(name, mass, density))
}

get_PT <- function(string){
  # Get Temperature
  start_pos <- gregexpr(text = string, pattern = "T = ")[[1]][1] + 4
  stop_pos <- gregexpr(text = string, pattern = "\\(C\\)")[[1]][1] - 2
  temperature <- substr(x = string, start = start_pos, stop = stop_pos)
  
  # Get Pressure
  start_pos <- gregexpr(text = string, pattern = "P = ")[[1]][1] + 4
  stop_pos <- gregexpr(text = string, pattern = "\\(kbars\\)")[[1]][1] - 2
  pressure <- substr(x = string, start = start_pos, stop = stop_pos)
  
  return(c(temperature, pressure))
}

print_status("Calculating volume parameters for the plot...")

ignore_phases <- c("Total solids", "System")
if(exclude_liquid){ignore_phases <- c(ignore_phases, "Liquid")}

melts <- read.delim(file = "melts.out", sep = c("\n"), stringsAsFactors = FALSE)[,1]

report <- data.frame(Name = "", Mass = "", Density = ""
                     , Temperature = "", Pressure = ""
                     , stringsAsFactors = FALSE
                     )

# Create blank variable to store state info
temp_i <- 0
press_i <- 0

# Scan through melts.out file and retrieve data ####
print_status("Parsing MELTS output file...")

for(i in 1:length(melts)){
  # Update state info
  if(grepl("T = ", melts[i])){
    new_PT <- get_PT(melts[i])
    temp_i <- new_PT[1]; press_i <- new_PT[2]
  }
  
  # Assign phase info
  if(grepl("  mass = ", melts[i])){
    print(paste("Scanning line ", i, " of ", length(melts), ".", sep = ""))
    report <- rbind(report, c(get_info(melts[i]), temp_i, press_i))
  }
}

# Get rid of dummy row
report <- report[-1,]

# Convert strings to numeric
report$Mass <- as.numeric(report$Mass)
report$Density <- as.numeric(report$Density)
report$Temperature <- as.numeric(report$Temperature)
report$Pressure <- as.numeric(report$Pressure)


# Calculate volume for each phase
Volume <- report$Mass/report$Density
report <- cbind(report, Volume)

# Assign steps, percent crystals, temperature, and pressure to plot
Step <- rep(0, nrow(report))
step_i <- 1

for(i in 1:nrow(report)){
  Step[i] <- step_i
  if(report$Name[i] == "System"){
    step_i <- step_i + 1
  }
}
report <- cbind(Step, report)

# Filter for "duplicate" phases ####

# I haven't updated this code to account for immiscibility
#   (e.g. two feldspars) because, at the initial time of
#   writing, immiscible phases were labeled identically.
#   Therefore, phases with the same name are combined into
#   the same polygon.

print_status("Filtering duplicate/immiscible phases...")

for(i in unique(report$Step)){ # going to do this for every time step
  step_i_rows <- which(report$Step == i) # Get the rows in the report variable for step_i
  
  for(j in 1:length(unique(report$Name[step_i_rows]))){
    name_j <- report$Name[step_i_rows[j]]
    
    # Skip over entries that aren't important phases
    if(name_j %in% ignore_phases){next}
    
    # Combine "duplicate" phases into first occurrence of phase
    flag_j <- which(report$Name[step_i_rows] == name_j)
    if(length(flag_j) > 2){stop("More than two phases with the same name!")}
    if(length(flag_j) > 1){
      print(paste("Step", i, ": Averaging", name_j, "properties.", sep = " "))
      # add masses together
      mass_1 <- report$Mass[step_i_rows][flag_j[1]]
      mass_2 <- report$Mass[step_i_rows][flag_j[2]]
      report$Mass[step_i_rows][flag_j[1]] <- mass_1 + mass_2
      report$Mass[step_i_rows][flag_j[2]] <- 0
      # average (weighted) densities
      mass_ratio <- mass_1/(mass_1 + mass_2)
      rho_1 <- report$Density[step_i_rows][flag_j[1]]
      rho_2 <- report$Density[step_i_rows][flag_j[2]]
      report$Density[step_i_rows][flag_j[1]] <- (rho_1*mass_ratio + rho_2*(1-mass_ratio))
      report$Density[step_i_rows][flag_j[2]] <- 0
      
      # Recalculate volumes
      report$Volume[step_i_rows][flag_j[1]] <- report$Mass[step_i_rows][flag_j[1]] / report$Density[step_i_rows][flag_j[1]]
      report$Volume[step_i_rows][flag_j[2]] <- 0
    }
  }
}
flag <- numeric(0)
for(i in 1:nrow(report)){ # Clean up empty rows
  if(report$Mass[i] == 0){flag <- c(flag, i)}
}
if(length(flag) > 0){report <- report[-flag,]}

# Make sure phases always appear in the same order of the report ####
print_status("Sorting phases so they always stack in the same order...")

if(TRUE){
  order_key <- cbind.data.frame(1:length(unique(report$Name)), unique(report$Name))
  # order_key <- cbind.data.frame(1:length(phases), phases)
  colnames(order_key) <- c("Rank", "Name")
  order_key$Rank <- as.numeric(order_key$Rank)
  order_key$Name <- as.character(order_key$Name)
  
  for(step_i in unique(report$Step)){
    
    #print(step_i)
    # Debug : step_i <- 1
    # Debug : step_i <- 31
    step_i_rows <- which(report$Step == step_i)
    
    chunk_i <- report[step_i_rows,]
    match_i <- match(order_key$Name, chunk_i$Name)
    match_i <- match_i[!is.na(match_i)]
    rank_i <- rank(match_i)
    report[step_i_rows,] <- report[step_i_rows[rank_i],]
  }
}

# Calculate polygon areas for plot ####
print_status("Calculating volume parameters for plot...")

y_bottom <- rep(0, nrow(report))
y_top <- rep(0, nrow(report))

# Max amount of crystallization
xl_max <- max(report$Mass[report$Name == "Total solids"])
Crystallized_frac <- rep(0, nrow(report))

for(i in unique(report$Step)){ # going to do this for every time step
  step_i_rows <- which(report$Step == i) # Get the rows in the report variable for step_i
  
  y_cumulative <- 0 # reset stack to bottom
  
  # Assign percent crystallization
  flag_i <- which((report$Step == i) & (report$Name == "Total solids"))
  if(length(flag_i == 1)){
      Crystallized_frac[step_i_rows] <- report$Mass[flag_i]/xl_max
     }
  
  for(j in 1:length(unique(report$Name[step_i_rows]))){ # going to do this for eeach phase
    name_j <- report$Name[step_i_rows[j]]
    
    # Skip over entries that aren't important phases
    if(name_j %in% ignore_phases){next}
    
    # Set bottom of jth polygon on top of previous polygons...
    y_bottom[step_i_rows[j]] <- y_cumulative
    # ...and set the top to the bottom plus volume of jth phase
    y_top_j <- y_cumulative + report$Volume[step_i_rows[j]]
    y_top[step_i_rows[j]] <- y_top_j
    
    # reset cumulative volume to top of the phase
    y_cumulative <- y_top_j
  }
  
  if(norm_volume & (y_cumulative > 0)){
    k_norm <- 100/y_cumulative
    y_bottom[step_i_rows] <- k_norm * y_bottom[step_i_rows]
    y_top[step_i_rows] <- k_norm * y_top[step_i_rows]
  }
  
}

# Add plotting parameters to report dataframe
report <- cbind(report, Crystallized_frac, y_bottom, y_top)

# Draw plot ####

# Saves a .pdf version of the plot to working directory

print_status("Saving phase volume figure...")
pdf(file = "figure.pdf", width = 8, height = 6, useDingbats = FALSE)

plot(0, 0, type = "n"
     , xlim = c(min(report[,settings[3,2]]), max(report[,settings[3,2]])), xlab = settings[3,2]
     , ylim = c(0, max(y_top)), ylab = "Volume, cm^3"
     )

phases <- unique(report$Name)[which(!(unique(report$Name) %in% ignore_phases))]
phase_cols <- hsv(1:length(phases)/length(phases), 1, 0.7)

for(i in 1:length(phases)){
  flag_i <- which(report$Name == phases[i])
  
  xs_i <- report[flag_i, settings[3, 2]]
  xs_i <- c(xs_i, rev(xs_i))
  
  ys_i <- c(report$y_bottom[flag_i], rev(report$y_top[flag_i]))
  
  polygon(xs_i, ys_i, col = phase_cols[i], lwd = 0.5)
  text(mean(xs_i), mean(ys_i), adj = c(0.5, 0.5)
       , labels = phases[i]
       )
}

dev.off()

# Plot key ####

# This section is probably unnecessary after adding labels to the plot!

print_status("Saving color key...")
pdf(file = "key.pdf", width = 4, height = 4, useDingbats = FALSE)
plot(0, 0, type = "n"
     , xlim = c(0, 10), ylim = c(0, length(phases)+1)
     , xlab = "", ylab = "", axes = FALSE
     )
for(i in 1:length(phases)){
  text(5, i, labels = phases[i], adj = c(0.5, 0.5)
       , col = phase_cols[i], cex = 1.5, lwd = 2
       )
}
dev.off()

# Save summary report ####

# Table version of data used in the plots

print_status("Saving summary report...")
write.csv(x = report, file = "report.csv", quote = FALSE, row.names = FALSE)

if(make_beep){beep(1)}
print_status("All done!")
