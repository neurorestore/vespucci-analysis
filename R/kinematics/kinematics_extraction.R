library(tidyverse)
library(magrittr)

## Load data
find_files_helper_function = function(input_dir, dir_type, pattern=NA){
	files = list.files(input_dir, full.names=T)
	
	## set left and right pattern to look for
	left_pattern = ifelse(dir_type == 'gait', '_L.csv', '_LDLC.*')
	right_pattern = ifelse(dir_type == 'gait', '_R.csv', '_RDLC.*')
	
	left_colname = paste0('left_', dir_type, '_file')
	right_colname = paste0('right_', dir_type, '_file')
	
	## match left and right files
	left_files = data.frame(
		left_file = files[grepl(left_pattern, files)] ) %>%
		mutate(id = gsub(left_pattern, '', basename(left_file))) %>%
		set_colnames(c(left_colname, 'id'))
	right_files = data.frame(
		right_file=files[grepl(right_pattern, files)]
		) %>%
		mutate(id = gsub(right_pattern, '', basename(right_file)))%>%
		set_colnames(c(right_colname, 'id'))
		
	files_df = 
		full_join(left_files, right_files, by='id') %>%
		relocate(id, .before=!!left_colname)
	files_df
}

## Create different gait event from DLC and Gait files	
full_wrapper = function(gait_file, dlc_file, foot, continuous_gait_events=T){
	print(gait_file)
	print(dlc_file)
	print(foot)
	## read gait file
	gait_df = read.csv(gait_file) %>%
		dplyr::select(Slice, Counter) %>%
		set_colnames(c('frame_index', 'event_code')) %>%
		mutate(
			foot = !!foot,
			event_name = as.character(event_code),
			event_name =fct_recode(
				event_name,
				'foot_strike' = '0',
				'drag' = '1',
				'foot_off' = '2'
			)
		) %>%
		distinct()
	## read dlc file
	dlc_header = read.csv(dlc_file, nrows=2)
	dlc_header = paste0(dlc_header[1,], '_', dlc_header[2,])
	dlc_df = read.csv(dlc_file, skip = 3) %>%
		set_colnames(dlc_header) %>%
		dplyr::rename_all(tolower) %>%
		dplyr::rename(frame_idx = 1)
	if("angle_x" %in% colnames(dlc_df) || "angle_y" %in% colnames(dlc_df)){
		colnames(dlc_df)[which(colnames(dlc_df) == "angle_x")] = "ankle_x"
		colnames(dlc_df)[which(colnames(dlc_df) == "angle_y")] = "ankle_y"
	}
	
	## creat pseudo event if the gait file is empty (animal do not walk)
	if (nrow(gait_df) == 1){
		gait_df = CreatPseudoEvents(foot, dlc_df, number_of_events = 10, frames_to_trim=10)
		gait_df = gait_df %>%
			dplyr::select(Slice, Counter) %>%
			set_colnames(c('frame_index', 'event_code')) %>%
			mutate(
				foot = !!foot,
				event_name = as.character(event_code),
				event_name =fct_recode(
					event_name,
					'foot_strike' = '0',
					'drag' = '1',
					'foot_off' = '2'
				)
			)
	}
	
	## get events
	check = 0
	event_idx = 1
	new_event = T
	event_finish = F
	gait_event_list = list()
	gait_event_error_log = data.frame()
	prev_event_name = 'foot_strike'
	prev_dlc_idx = gait_df$frame_index[1]
	first_event = TRUE

	for (i in 1:nrow(gait_df)){
		event_code = gait_df$event_code[i]
		# print(i)
		# print(new_event)
		# print(check)
		# print(event_code)
		if (event_code == 0){
			new_event = T
			start_idx = i
			start_frame_idx = gait_df$frame_index[i]
			end_idx = NULL
		}
		
		# fix for mixed 0-0 and 0-1
		if (!continuous_gait_events & event_code == 1 & check == 0 & !first_event & new_event){
				if (!is.null(end_idx)) {
						start_idx = end_idx
						start_frame_idx = end_frame_idx
				}
				check = 1
				# print('Here')
		}
		
		## check if cycle is completed
		if (
			(event_code != check & new_event) ||
			(is.na(gait_df$event_code[i+3]) & new_event & gait_df$event_code[i] == 0) || 
			(any(duplicated((gait_df$frame_index[i:(i+3)]))) & !(is.na(gait_df$event_code[i+3])))
		) {
			print(i)
			if (event_code != check & new_event){
					error_msg = paste0('Expecting ', check, ' but got ', event_code)
			} else if ((is.na(gait_df$event_code[i+3]) & new_event & gait_df$event_code[i] == 0)) {
					error_msg = 'Error 2'
			} else if (any(duplicated((gait_df$frame_index[i:(i+3)]))) & !(is.na(gait_df$event_code[i+3]))) {
					error_msg = 'Error 3'
			}
			print(error_msg)
			error = data.frame(
				'foot'=foot,
				'event_idx'=event_idx,
				'index'=i, 
				'error'=error_msg
			)
			gait_event_error_log = rbind(gait_event_error_log, error)
			check = 0
			new_event = F
			event_idx =  event_idx + 1
		}
		
		if (new_event){
			check = check + 1
			if (event_finish & event_code == 0) {
				check = 0
				event_finish = F
			}

			if (!continuous_gait_events & check > 2) {
				event_finish = T
			}

			if (check > 2) {
				check = 0
				end_idx = i+1
				end_frame_idx = gait_df$frame_index[i+1]
				
				event_df = gait_df[start_idx:end_idx,] %>%
					mutate(event_idx = event_idx) %>%
					set_rownames(NULL)
				
				event_coords_df = data.frame()
				
				if (!any(is.na(event_df))){
						for (j in 1:(nrow(event_df))){
								curr_frame_idx = event_df$frame_index[j]
								curr_dlc_idx = which(dlc_df$frame_idx == curr_frame_idx)
								curr_event_name = event_df$event_name[j]
								
								if (prev_dlc_idx != curr_dlc_idx & prev_dlc_idx != gait_df$frame_index[nrow(gait_df)]){
										event_coords_name = paste0(prev_event_name, '_to_', curr_event_name)
										if(j == nrow(event_df)){
												temp_event_coords_df = dlc_df[prev_dlc_idx:curr_dlc_idx,] %>%
														set_rownames(NULL) %>%
														mutate(event_name = event_coords_name)
												event_coords_df = rbind(event_coords_df, temp_event_coords_df)
										} else{
												temp_event_coords_df = dlc_df[prev_dlc_idx:(curr_dlc_idx-1),] %>%
														set_rownames(NULL) %>%
														mutate(event_name = event_coords_name)
												event_coords_df = rbind(event_coords_df, temp_event_coords_df)
										}
										
								}
								
								prev_frame_idx = curr_frame_idx
								prev_dlc_idx = curr_dlc_idx
								prev_event_name = curr_event_name
						}
						
						
						event_list = list(
								'error_log' = data.frame(
										'Check'='Event cycle',
										'Error'='None'
								),
								'event_df'=event_df,
								'event_coords_df'=event_coords_df
						)
						
						gait_event_list[[paste0('Event_', event_idx)]] = event_list
						event_idx = event_idx + 1   
						first_event = FALSE
						}
				}
		}
	}
	
	gait_event_list
}

## Check event likelihood
check_gait_event_likelihood = function(event, threshold) {
	
	likelihood_cols = colnames(event$event_coords_df)
	likelihood_cols = likelihood_cols[grepl('likelihood', likelihood_cols)]
	
	check = unlist(sapply(1:nrow(event$event_df), function(x){
		
		frame_idx = event$event_df$frame_index[x]
		event_name = event$event_df$event_name[x]
		
		frame_likelihoods = event$event_coords_df %>%
			filter(frame_idx == !!frame_idx) %>%
			dplyr::select(likelihood_cols)
		frame_likelihoods = names(frame_likelihoods)[frame_likelihoods < threshold]
		if (length(frame_likelihoods) > 0) {
			frame_likelihoods = paste0(event_name, '--', frame_likelihoods)
		} 
		frame_likelihoods
	}))
	check = paste0(check, collapse = '; ')
	error = 'None'
	if (check != ''){
		error = check
	}
	event$error_log = rbind(
		event$error_log,
		data.frame('Check'='Event frame failed qc', 'Error'=error))
	event
}

## Filter event with an high likelihood
filter_event_likehood = function(event, threshold, minimum_event_frames_percentage) {
	
	event_coords_df = event$event_coords_df
	event_df = event$event_df
	
	likelihood_cols = colnames(event_coords_df)
	likelihood_cols = likelihood_cols[grepl('likelihood', likelihood_cols)]
	
	filter = rowSums(event_coords_df[,likelihood_cols] < threshold) == 0
	filtered_event_coords_df = event_coords_df[filter,] %>%
		set_rownames(NULL) %>%
		dplyr::select(-likelihood_cols)
	
	event$filtered_event_coords_df = filtered_event_coords_df
	
	error = 'None'
	if (mean(filter) < minimum_event_frames_percentage){
		error = 'Not enough frames'
	}
	
	event$error_log = rbind(
		event$error_log,
		data.frame('Check'='Number of frames after qc', 'Error'=error)
	)
	
	event
}

## Reverse DLC data
reverse_data = function(event) {
	event$filtered_event_coords_df$crest_x = -event$filtered_event_coords_df$crest_x
	event$filtered_event_coords_df$crest_y = -event$filtered_event_coords_df$crest_y
	event$filtered_event_coords_df$hip_x = -event$filtered_event_coords_df$hip_x
	event$filtered_event_coords_df$hip_y = -event$filtered_event_coords_df$hip_y
	event$filtered_event_coords_df$knee_x = -event$filtered_event_coords_df$knee_x
	event$filtered_event_coords_df$knee_y = -event$filtered_event_coords_df$knee_y
	event$filtered_event_coords_df$ankle_x = -event$filtered_event_coords_df$ankle_x
	event$filtered_event_coords_df$ankle_x = -event$filtered_event_coords_df$ankle_x
	event$filtered_event_coords_df$foot_x = -event$filtered_event_coords_df$foot_x
	event$filtered_event_coords_df$foot_y = -event$filtered_event_coords_df$foot_y
	
	event
}

## Estimate the real y position
estimate_runway_y_position = function(event, x_position){
	features = event$features
	return(
		(features$estimated_runway_slope * x_position) + features$estimated_runway_constant)
}

## Get main features needed for kinematics ectaction
get_meta_features_function = function(event, frame_rate, pixel_value) {
	event_coords_df = event$filtered_event_coords_df
	features = list()
	start_idx = 1
	end_idx = nrow(event_coords_df)
	
	features$start_idx = start_idx
	features$end_idx = end_idx
	features$start_time = event_coords_df$frame_idx[start_idx] / frame_rate
	features$end_time = event_coords_df$frame_idx[end_idx] / frame_rate
	features$estimated_runway_slope = 
		(event_coords_df$foot_y[end_idx] - event_coords_df$foot_y[start_idx])/
		(event_coords_df$foot_x[end_idx] - event_coords_df$foot_x[start_idx])
	features$estimated_runway_constant =  
		event_coords_df$foot_y[end_idx] -
		(features$estimated_runway_slope * event_coords_df$foot_x[end_idx])
	features$frame_rate = frame_rate
	features$pixel_value = pixel_value
	
	event$features = features
	event
}

calculate_2d_euclidean_distance = function(start_coords, end_coords) {
	dist = sqrt(((start_coords[1] - end_coords[1])**2)+((start_coords[2] - end_coords[2])**2))
	return (dist)
}

## Calculate the different needed angles
calculate_hip_angle = function(event){
	hip_crest_x = event$filtered_event_coords_df$crest_x - event$filtered_event_coords_df$hip_x
	hip_crest_y = ((event$filtered_event_coords_df$crest_y) - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x)) -
		((event$filtered_event_coords_df$hip_y) - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))
	knee_hip_x = event$filtered_event_coords_df$knee_x - event$filtered_event_coords_df$hip_x
	knee_hip_y = ((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
		((event$filtered_event_coords_df$hip_y) - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))
	
	dot_prod = c()
	for (i in 1:length(hip_crest_x)){
		a = c(hip_crest_x[i], hip_crest_y[i])
		b = c(knee_hip_x[i], knee_hip_y[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_knee_angle = function(event){
	knee_hip_x = event$filtered_event_coords_df$hip_x - event$filtered_event_coords_df$knee_x
	knee_hip_y = ((event$filtered_event_coords_df$hip_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$hip_x)) -
		((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$knee_x))
	ankle_knee_x = event$filtered_event_coords_df$ankle_x - event$filtered_event_coords_df$knee_x
	ankle_knee_y = ((event$filtered_event_coords_df$ankle_x) - estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x)) -
		((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$knee_x))
	
	dot_prod = c()
	for (i in 1:length(knee_hip_x)){
		a = c(knee_hip_x[i], knee_hip_y[i])
		b = c(ankle_knee_x[i], ankle_knee_y[i])
		angle = (atan2(b[2],b[1]) - atan2(a[2], a[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_ankle_angle = function(event){
	ankle_knee_x = event$filtered_event_coords_df$knee_x - event$filtered_event_coords_df$ankle_x
	ankle_knee_y = ((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$knee_x))-
		((event$filtered_event_coords_df$ankle_x) - estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x))
	foot_ankle_x = event$filtered_event_coords_df$foot_x - event$filtered_event_coords_df$ankle_x
	foot_ankle_x = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))-
		((event$filtered_event_coords_df$ankle_x) - estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x))
	
	dot_prod = c()
	for (i in 1:length(ankle_knee_x)){
		a = c(ankle_knee_x[i], ankle_knee_y[i])
		b = c(foot_ankle_x[i], foot_ankle_x[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_limb_angle = function(event){
	crest_ground_x = event$filtered_event_coords_df$crest_x - event$filtered_event_coords_df$crest_x
	crest_ground_y = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))-
		((event$filtered_event_coords_df$crest_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$crest_x))
	foot_crest_x = event$filtered_event_coords_df$foot_x - event$filtered_event_coords_df$crest_x
	foot_crest_y = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x)) -
		((event$filtered_event_coords_df$crest_y)- estimate_runway_y_position(event, event$filtered_event_coords_df$crest_x))
	dot_prod = c()
	for (i in 1:length(foot_crest_x)){
		a = c(crest_ground_x[i], crest_ground_y[i])
		b = c(foot_crest_x[i], foot_crest_y[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_MTP_angle = function(event){
	foot_ground_x = event$filtered_event_coords_df$foot_x - event$filtered_event_coords_df$foot_x
	foot_ground_y = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x) - 10)-
		((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))
	ankle_foot_x = event$filtered_event_coords_df$ankle_x - event$filtered_event_coords_df$foot_x
	ankle_foot_y = ((event$filtered_event_coords_df$ankle_x) - estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x))-
		((event$filtered_event_coords_df$foot_y)-estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))
	dot_prod = c()
	for (i in 1:length(ankle_foot_x)){
		a = c(foot_ground_x[i], foot_ground_y[i])
		b = c(ankle_foot_x[i], ankle_foot_y[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_thigh_angle = function(event){
	knee_hip_x = event$filtered_event_coords_df$knee_x - event$filtered_event_coords_df$hip_x
	knee_hip_y = ((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$knee_x)) -
		((event$filtered_event_coords_df$hip_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$hip_x))
	ground_hip_x = event$filtered_event_coords_df$hip_x - event$filtered_event_coords_df$hip_x
	ground_hip_y = ((event$filtered_event_coords_df$knee_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$knee_x)) -
		((event$filtered_event_coords_df$hip_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$hip_x))
	dot_prod = c()
	for (i in 1:length(ground_hip_x)){
		a = c(knee_hip_x[i], knee_hip_y[i])
		b = c(ground_hip_x[i], ground_hip_y[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_crest_angle = function(event){
	crest_hip_x = event$filtered_event_coords_df$crest_x - event$filtered_event_coords_df$hip_x
	crest_hip_y = ((event$filtered_event_coords_df$crest_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$crest_x))-
		((event$filtered_event_coords_df$hip_y)-estimate_runway_y_position(event, event$filtered_event_coords_df$hip_x))
	ground_hip_x = event$filtered_event_coords_df$hip_x - event$filtered_event_coords_df$hip_x
	ground_hip_y = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))-
		((event$filtered_event_coords_df$hip_y)-estimate_runway_y_position(event, event$filtered_event_coords_df$hip_x))
	dot_prod = c()
	for (i in 1:length(ground_hip_x)){
		a = c(crest_hip_x[i], crest_hip_y[i])
		b = c(ground_hip_x[i], ground_hip_y[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_shank_angle = function(event){
	ankle_knee_x = event$filtered_event_coords_df$ankle_x - event$filtered_event_coords_df$knee_x
	ankle_knee_y = ((event$filtered_event_coords_df$ankle_x)-estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x))-
		((event$filtered_event_coords_df$knee_y)-estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x))
	ground_knee_x = event$filtered_event_coords_df$knee_x - event$filtered_event_coords_df$knee_x
	ground_knee_y = ((event$filtered_event_coords_df$foot_y)-estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))-
		((event$filtered_event_coords_df$knee_y)-estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x))
	
	dot_prod = c()
	for (i in 1:length(ankle_knee_x)){
		a = c(ankle_knee_x[i], ankle_knee_y[i])
		b = c(ground_knee_x[i], ground_knee_y[i])
		angle = (atan2(b[2],b[1]) - atan2(a[2], a[1]))*( 180.0 / pi )
		if (angle < 0){
			angle = angle + 360
		}
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

calculate_foot_angle = function(event){
	foot_ankle_x = event$filtered_event_coords_df$foot_x - event$filtered_event_coords_df$ankle_x
	foot_ankle_x = ((event$filtered_event_coords_df$foot_y) - estimate_runway_y_position(event, event$filtered_event_coords_df$foot_x))-
		((event$filtered_event_coords_df$ankle_x)-event$filtered_event_coords_df$ankle_x)
	ground_ankle_x = event$filtered_event_coords_df$ankle_x - event$filtered_event_coords_df$ankle_x
	ground_ankle_x = ((event$filtered_event_coords_df$ankle_x)-estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x) - 10)-
		((event$filtered_event_coords_df$ankle_x)-estimate_runway_y_position(event, event$filtered_event_coords_df$ankle_x))
	dot_prod = c()
	for (i in 1:length(ground_ankle_x)){
		a = c(foot_ankle_x[i], foot_ankle_x[i])
		b = c(ground_ankle_x[i], ground_ankle_x[i])
		angle = (atan2(a[2],a[1]) - atan2(b[2], b[1]))*( 180.0 / pi )
		dot_prod = c(dot_prod, angle)
	}
	return(dot_prod)
}

get_general_gait_features_function = function(event) {
	
	analysis = list()
		
	features = event$features
	frame_rate = features$frame_rate
	pixel_value = features$pixel_value
	start_idx = features$start_idx
	end_idx = features$end_idx
	event_coords_df = event$filtered_event_coords_df
	
	## get gait cycle basic features
	
	analysis$gait_cycle_onset = event_coords_df$frame_idx[start_idx]
	analysis$gait_cycle_end = event_coords_df$frame_idx[end_idx]
	analysis$gait_cycle_duration = (analysis$gait_cycle_end - analysis$gait_cycle_onset)/ frame_rate
	
	## get gait events indices
	
	analysis$stance_start_idx = event_coords_df$frame_idx[start_idx]
	analysis$stance_end_idx = event_coords_df %>%
		filter(event_name == 'foot_strike_to_drag') %>%
		pull(frame_idx) %>%
		max()
	
	analysis$drag_start_idx = event_coords_df %>%
		filter(event_name == 'drag_to_foot_off') %>%
		pull(frame_idx) %>%
		min()
	
	analysis$drag_end_idx = event_coords_df %>%
		filter(event_name == 'drag_to_foot_off') %>%
		pull(frame_idx) %>%
		max()
	
	analysis$swing_start_idx = event_coords_df %>%
		filter(event_name == 'foot_off_to_foot_strike') %>%
		pull(frame_idx) %>%
		min()
	
	analysis$swing_end_idx = event_coords_df %>%
		filter(event_name == 'foot_off_to_foot_strike') %>%
		pull(frame_idx) %>%
		max()

	
	## finally generate some gait features
	analysis$stance_duration = ((analysis$stance_end_idx - analysis$stance_start_idx) / frame_rate)
	analysis$swing_duration = ((analysis$swing_end_idx - analysis$swing_start_idx) / frame_rate)
	analysis$swing_duration_percentage = (analysis$swing_duration/analysis$gait_cycle_duration) * 100
	analysis$drag_duration = ((analysis$drag_end_idx - analysis$drag_start_idx) / frame_rate)
	analysis$drag_duration_percentage = (analysis$drag_duration/analysis$gait_cycle_duration) * 100
	
	
	analysis$stride_length = (abs(
		event_coords_df %>% 
			filter(frame_idx == analysis$swing_end_idx) %>%
			pull(ankle_x) %>% 
			unique() - 
			event_coords_df %>% 
			filter(frame_idx == analysis$stance_end_idx) %>%
			pull(ankle_x) %>% 
			unique()
	)/pixel_value)
	analysis$step_length = (abs(
		event_coords_df %>% 
			filter(frame_idx == analysis$swing_end_idx) %>%
			pull(ankle_x) %>% 
			unique() - 
		event_coords_df %>% 
			filter(frame_idx == analysis$swing_start_idx) %>%
			pull(ankle_x) %>% 
			unique()
	)/pixel_value)
	
	analysis$path_length = (calculate_2d_euclidean_distance(
		event_coords_df %>%
			filter(frame_idx == analysis$swing_start_idx) %>%
			dplyr::select(ankle_x, ankle_y) %>%
			distinct() %>%
			unlist(),
		event_coords_df %>%
			filter(frame_idx == analysis$swing_end_idx) %>%
			dplyr::select(ankle_x, ankle_y) %>%
			distinct() %>%
			unlist()
	)[[1]]/pixel_value)
	
	analysis$animal_size_factor =
		mean(sapply(1:nrow(event_coords_df), function(x){
			calculate_2d_euclidean_distance(
			unlist(event_coords_df[x, c('crest_x', 'crest_y')]),
			unlist(event_coords_df[x, c('hip_x', 'hip_y')])
		)}))
	
	analysis$step_height_foot = (max(event_coords_df$foot_y) - 
								event_coords_df %>%
									filter(frame_idx == analysis$stance_end_idx | frame_idx == analysis$swing_end_idx) %>%
									pull(foot_y) %>%
									min())/pixel_value
	analysis$step_height_foot_norm = analysis$step_height_foot / analysis$animal_size_factor
	analysis$step_height_ankle = ((max(event_coords_df$ankle_y) - 
								event_coords_df %>%
									filter(frame_idx == analysis$stance_end_idx | frame_idx == analysis$swing_end_idx) %>%
									pull(ankle_y) %>%
									min())/pixel_value)
	analysis$step_height_ankle_norm = analysis$step_height_ankle / analysis$animal_size_factor
	
	analysis$min_elev_crest = ((abs((event$filtered_event_coords_df$crest_y - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x)) -
																	 (event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))) %>%
		min())/pixel_value)
	analysis$max_elev_crest = ((abs((event$filtered_event_coords_df$crest_y - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x)) -
																	(event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))) %>%
		max())/pixel_value)
	
	analysis$min_elev_thigh = ((abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																	 (event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))) %>%
		min())/pixel_value)
	analysis$max_elev_thigh = ((abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																	(event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))) %>%
		max())/pixel_value)
	
	analysis$min_elev_shank = ((abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																	 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))) %>%
		min())/pixel_value)
	analysis$max_elev_shank = ((abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																	(event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))) %>%
		max())/pixel_value)
	
	
	analysis$min_elev_foot = ((abs((event$filtered_event_coords_df$foot_y - estimate_runway_y_position(event,event$filtered_event_coords_df$foot_x)) -
																	 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))) %>%
		min())/pixel_value)
	analysis$max_elev_foot = ((abs((event$filtered_event_coords_df$foot_y - estimate_runway_y_position(event,event$filtered_event_coords_df$foot_x)) -
																 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))) %>%
		max())/pixel_value)
	
	
	analysis$min_elev_limb = ((abs((event$filtered_event_coords_df$foot_y - estimate_runway_y_position(event,event$filtered_event_coords_df$foot_x)) -
																	(event$filtered_event_coords_df$crest_y - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x))) %>%
		min())/pixel_value)
	analysis$max_elev_limb = ((abs((event$filtered_event_coords_df$foot_y - estimate_runway_y_position(event,event$filtered_event_coords_df$foot_x)) -
																 (event$filtered_event_coords_df$crest_y - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x))) %>%
		max())/pixel_value)
		
	analysis$ampl_elev_crest = analysis$max_elev_crest - analysis$min_elev_crest
	analysis$ampl_elev_thigh = analysis$max_elev_thigh - analysis$min_elev_thigh
	analysis$ampl_elev_shank = analysis$max_elev_shank - analysis$min_elev_shank
	analysis$ampl_elev_foot = analysis$max_elev_foot - analysis$min_elev_foot
	analysis$ampl_elev_limb = analysis$max_elev_foot - analysis$min_elev_limb
	
	analysis$max_hip_angle = max(calculate_hip_angle(event))
	analysis$min_hip_angle = min(calculate_hip_angle(event))
	analysis$ampl_hip_angle = analysis$max_hip_angle - analysis$min_hip_angle
	analysis$max_knee_angle = max(calculate_knee_angle(event))
	analysis$min_knee_angle = min(calculate_knee_angle(event))
	analysis$ampl_knee_angle = analysis$max_knee_angle - analysis$min_knee_angle
	analysis$max_ankle_angle = max(calculate_ankle_angle(event))
	analysis$min_ankle_angle = min(calculate_ankle_angle(event))
	analysis$ampl_ankle_angle = analysis$max_ankle_angle - analysis$min_ankle_angle
	
	analysis$min_angle_speed_hip = min(diff(calculate_hip_angle(event))/(1/event$features$frame_rate))
	analysis$min_angle_speed_knee = min(diff(calculate_knee_angle(event))/(1/event$features$frame_rate))
	analysis$min_angle_speed_ankle = min(diff(calculate_ankle_angle(event))/(1/event$features$frame_rate))
	analysis$min_angle_speed_limb = min(diff(calculate_limb_angle(event))/(1/event$features$frame_rate))
	analysis$max_angle_speed_hip = max(diff(calculate_hip_angle(event))/(1/event$features$frame_rate))
	analysis$max_angle_speed_knee = max(diff(calculate_knee_angle(event))/(1/event$features$frame_rate))
	analysis$max_angle_speed_ankle = max(diff(calculate_ankle_angle(event))/(1/event$features$frame_rate))
	analysis$max_angle_speed_limb = max(diff(calculate_limb_angle(event))/(1/event$features$frame_rate))
	analysis$ampl_angle_speed_hip = analysis$max_angle_speed_hip - analysis$min_angle_speed_hip
	analysis$ampl_angle_speed_knee = analysis$max_angle_speed_knee  - analysis$min_angle_speed_knee
	analysis$ampl_angle_speed_ankle = analysis$max_angle_speed_ankle - analysis$min_angle_speed_ankle
	analysis$ampl_angle_speed_limb = analysis$max_angle_speed_limb - analysis$min_angle_speed_limb
	
	analysis$corr_crest_thigh = cor(abs((event$filtered_event_coords_df$crest_y - estimate_runway_y_position(event,event$filtered_event_coords_df$crest_x)) -
																					(event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))),
																		abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																					(event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))))
	
	analysis$corr_thight_shank = cor(abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																				 (event$filtered_event_coords_df$hip_y - estimate_runway_y_position(event,event$filtered_event_coords_df$hip_x))),
																	 abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																				 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))))
	
	analysis$corr_shank_foot = cor(abs((event$filtered_event_coords_df$knee_y - estimate_runway_y_position(event,event$filtered_event_coords_df$knee_x)) -
																			 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))),
																 abs((event$filtered_event_coords_df$foot_y - estimate_runway_y_position(event,event$filtered_event_coords_df$foot_x)) -
																			 (event$filtered_event_coords_df$ankle_x - estimate_runway_y_position(event,event$filtered_event_coords_df$ankle_x))))
	
	analysis$corr_hip_knee_joints = cor(calculate_hip_angle(event), calculate_knee_angle(event))
	analysis$corr_knee_ankle_joints = cor(calculate_knee_angle(event), calculate_ankle_angle(event))
	analysis$corr_ankle_foot_joints = cor(calculate_ankle_angle(event), calculate_MTP_angle(event))
	
	analysis$dur_pos_bw_crest_thigh = (which.max(calculate_crest_angle(event)) - which.min(calculate_thigh_angle(event)))*100/event$features$end_idx
	analysis$dur_pos_fw_crest_thigh = (which.min(calculate_crest_angle(event)) - which.max(calculate_thigh_angle(event)))*100/event$features$end_idx
	analysis$dur_pos_bw_thigh_shank = (which.min(calculate_thigh_angle(event)) - which.max(calculate_shank_angle(event)))*100/event$features$end_idx
	analysis$dur_pos_fw_thigh_shank = (which.max(calculate_thigh_angle(event)) - which.min(calculate_shank_angle(event)))*100/event$features$end_idx
	analysis$dur_pos_bw_shank_foot = (which.max(calculate_shank_angle(event)) - which.min(calculate_foot_angle(event)))*100/event$features$end_idx
	analysis$dur_pos_fw_shank_foot = (which.min(calculate_shank_angle(event)) - which.max(calculate_foot_angle(event)))*100/event$features$end_idx
	
	crest_hip_diff = event$event_coords_df$crest_x - event$event_coords_df$ankle_x
	
	analysis$max_bw_pos = max(crest_hip_diff)/pixel_value
	analysis$max_fw_pos = min(crest_hip_diff)/pixel_value
	analysis$bw_fw_amp = abs(analysis$max_bw_pos - analysis$max_fw_pos)
	
	max_y_hip = event$event_coords_df %>%
		arrange(-hip_y) %>%
		dplyr::select(hip_y, foot_x, frame_idx) %>%
		dplyr::slice(1)
	analysis$max_vertical_hip_displacement = 
		max_y_hip$hip_y - estimate_runway_y_position(event, max_y_hip$foot_x)
	
	event$analysis = list(
		'general_gait_features' = analysis
	)
	return (event)
}

## Create pseudo event for non-walking animal
CreatPseudoEvents = function(foot, dlc_df, number_of_events, frames_to_trim=10){
	output_df = data.frame()
	event_frame_range = c()
	if(foot == "left"){
		sequence_start = frames_to_trim
		original_sequence_end = nrow(dlc_df)-frames_to_trim
		adjusted_sequence_end = original_sequence_end - ((original_sequence_end - sequence_start + 1) %% 10)
		sequence_to_split = seq(sequence_start,adjusted_sequence_end )
		subarray_length = length(sequence_to_split) %/% number_of_events
		split_factor = gl(number_of_events, subarray_length)
		event_frame_range = split(sequence_to_split, split_factor)
	} else if (foot == "right"){
		sequence_start = frames_to_trim
		original_sequence_end = nrow(dlc_df)-frames_to_trim
		adjusted_sequence_end = original_sequence_end - ((original_sequence_end - sequence_start + 1) %% 10)
		sequence_to_split = seq(sequence_start,adjusted_sequence_end)
		subarray_length = length(sequence_to_split) %/% number_of_events
		split_factor = gl(number_of_events, subarray_length)
		event_frame_range = split(sequence_to_split, split_factor)
	}
	for (i in 1:length(event_frame_range)){
		curr_event_frame_range = event_frame_range[i]
		start_frame = curr_event_frame_range[[1]][1]
		end_frame = curr_event_frame_range[[1]][length(curr_event_frame_range[[1]])]
		
		temp_df = data.frame(
			Slice = c(start_frame, start_frame+1, end_frame),
			Counter = c(0, 1, 2)
		)
		output_df = rbind(output_df, temp_df)
	}
	rownames(output_df) = NULL
	return(output_df)
}

## Extract parameters for a single species #######################
extract_parameters_1species_function = function(species_path, feet, conditions, species_name, frame_rate, pixel_to_cm, indiv_name = "no_name", continuous_gait_events=F){
	
	## Initialize the final structure of the master file
	## BE CAREFUL: change ncol if you add or supress some analysis parameters
	store_features = data.frame()
	for(j in 1:length(conditions)){
		# set the directory for the different studied conditions
		condition_dir = paste(species_path, "/" ,conditions[j], sep="")
		print(condition_dir)
		## get matching gait files
		gait_files_df = find_files_helper_function(condition_dir, 'gait')
		## get matching dlc files
		dlc_files_df = find_files_helper_function(condition_dir, 'DLC')
		## find matching gait + dlc files
		files_df = full_join(gait_files_df, dlc_files_df)
		
		## get missing files
		rows_with_missing_files = rowSums(is.na(files_df)) > 0
		missing_files_log = files_df[rows_with_missing_files,]
		files_df = files_df[!rows_with_missing_files,]
		
		## finally start analysing
		threshold = 0.9
		frame_rate = frame_rate
		pixel_to_cm = pixel_to_cm
		minimum_event_frames_percentage = 0.8
		
		## iterate on the different files of a given condition
		for(k in 1:nrow(files_df)){
			file_args = files_df[k,]
			print(file_args)
			
			## iterate on the two feet
			for(l in 1:length(feet)){
				foot_gait_file = paste(feet[l], "_gait_file", sep="")
				foot_dlc_file = paste(feet[l], "_DLC_file", sep="")
				gait_file = file_args[[foot_gait_file]]
				dlc_file = file_args[[foot_dlc_file]]
				
				gait_event_list = full_wrapper(gait_file, dlc_file, feet[l], continuous_gait_events=T)
				
				for(i in 1:length(gait_event_list)){
					event = gait_event_list[[i]]
					event = check_gait_event_likelihood(event, threshold = threshold)
					event$error_log
					event = filter_event_likehood(event, threshold, minimum_event_frames_percentage)
					if (all(event$error_log$Error=='None')){
						event = reverse_data(event)
						event = get_meta_features_function(event, frame_rate, pixel_to_cm)
						event = get_general_gait_features_function(event)
						replicate = k
						if(indiv_name != "no_name"){
							replicate = indiv_name
						}
						row = paste(species_name, conditions[j], replicate ,feet[l], 'event', i)
						event$analysis$general_gait_features = lapply(event$analysis$general_gait_features, function(x) ifelse(length(x) == 0, 0, x))
						df = data.frame(event$analysis$general_gait_features)
						colnames(df) = colnames(store_features)
						rownames(df) = row
						store_features = rbind(store_features, df)
						colnames(store_features) = names(event$analysis$general_gait_features)
					}
				}
			}
		}
	}
	return(store_features)
}

## Extract parameters for all species 
extract_parameters_function = function(data_path, output_dir, config_file, conditions, continuous_gait_events=T){
	
	read_config = read.csv(config_file)
	
	file_species_list = list.dirs(data_path, recursive = F)
	
	## Initialize the final structure of the master file
	## BE CAREFUL: change ncol if you add or supress some analysis parameters
	store_features = data.frame(matrix(ncol = 74, nrow = 1))
	for(species_path in file_species_list){
		species = basename(species_path)
		frame_rate = read_config$Frame_rate[read_config$Species == species]
		pixel_to_cm = read_config$Pixel_value[read_config$Species == species]
		feet = c("right", "left")
		file_list = list.files(species_path)
		
		if (all(conditions %in% file_list)) {
			df = extract_parameters_1species_function(species_path, feet, conditions, species, frame_rate, pixel_to_cm, continuous_gait_events=T)
			colnames(store_features) = colnames(df)
			store_features = rbind(store_features, df)
		} else if (!(any(conditions %in% file_list))){
			list_indiv = list.files(getwd())
			for (indiv in 1:length(list_indiv)){
				indiv_path = paste("./", list_indiv[indiv], sep="")
				setwd(indiv_path)
				df = extract_parameters_1species_function(feet, conditions, file_species_list[species], frame_rate, pixel_to_cm, continuous_gait_events, indiv_name = list_indiv[indiv])
				colnames(store_features) = colnames(df)
				store_features = rbind(store_features, df)
				setwd("../")
			}
		} else{
			print("some conditions are missing")
		}
	 
	}
	store_features <- store_features[rowSums(is.na(store_features)) != ncol(store_features), ]
	write.csv(store_features, file = paste0(output_dir, "master_file.csv"))
	return(store_features)
}

#############################################################################
############ Call the function to create master file ########################
data_path = "data/kinematics/raw_data"
config_file = "data/kinematics/kinematics_config.csv"
output_dir = 'data/kinematices/'
conditions = basename(list.dirs(paste0(data_path, 'mouse/'), recursive=F))
extract_parameters_function(data_path, output_dir, config_file, conditions)

dat = read.csv("data/kinematics/master_file.csv", row.name=1) %>%
	mutate(file = row.names(.)) %>%
	separate(file, c("species", "condition", "replicate", "foot", "event", "event_no"), " ") %>%
	relocate("species", "condition", "replicate", "foot") %>%
	dplyr::select(-event, -event_no, -gait_cycle_onset, -gait_cycle_end, -stance_start_idx, -stance_end_idx, -drag_start_idx, -drag_end_idx, -swing_start_idx, -swing_end_idx) %>%
		drop_na()
dat %>% dplyr::select(condition, replicate) %>% table()
dat %<>% mutate(replicate = paste0(condition, '_', replicate))
saveRDS(dat, 'data/kinematics/final_data.rds')