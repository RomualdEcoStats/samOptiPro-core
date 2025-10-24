# Construction on appropriate inits
# SIMULATION from the nimble model with some nodes fixed to specified values
# 14/12/2021
# R. Patin

simul_inits <- function(Fix, Data_nimble, Const_nimble, modelfile,
												logfile_connection){
	
	# Load the model and create and object named mymod to be used for simulation ----------------
	hindcast <- TRUE
	source(modelfile)
	inits <- list()
	n.chains <- length(Fix)
	for(thischain in 1:n.chains){
		if(thischain == 1) {
			cli_process_start("Building model")
			sink(file = logfile_connection, append = TRUE)
			sink(file = logfile_connection,type = "message",append = TRUE)
			mymod <- nimbleModel(code = model.nimble, name = 'model.nimble', 
													 constants = Const_nimble, data = Data_nimble)
			sink()
			sink(type = "message")
			cli_process_done()
			
		}
		cli_process_start("Simulating inits for chain {thischain}")
		sink(file = logfile_connection, append = TRUE)
		sink(file = logfile_connection,type = "message",append = TRUE)
		
		# Assign fixed values (carefully chosen) to model nodes that shouldn't be simulated ---------
		
		
		mymod$N1_pr          <- Fix[[thischain]]$N1_pr
		mymod$N2_pr          <- Fix[[thischain]]$N2_pr
		mymod$N2          <- Fix[[thischain]]$N2
		mymod$log_N9_pr          <- log(Fix[[thischain]]$log_N9_pr)
		mymod$h6_hw       <- Fix[[thischain]]$h6_hw
		mymod$h9_hw       <- Fix[[thischain]]$h9_hw
		mymod$h5       <- Fix[[thischain]]$h5
		mymod$h8   <- Fix[[thischain]]$h8
		mymod$logit_theta3   <- Fix[[thischain]]$logit_theta3
		mymod$logit_theta4   <- Fix[[thischain]]$logit_theta4
		mymod$sd_theta3   <- Fix[[thischain]]$sd_theta3
		mymod$sd_theta4   <- Fix[[thischain]]$sd_theta4
		
		
		# Simulate all nodes of the model except those that have been fixed  -----------------
		# Step 1 - Identify the indices of nodes that must be simulated ---------------------
		#         (versus nodes that are fixed)	 
		# 
		
		# All Nodes of the model
		
		All_Nodes <- mymod$getNodeNames(returnScalarComponents = TRUE)
		head(All_Nodes)
		All_Nodes_names <- mymod$getVarNames(includeLogProb = FALSE, nodes=All_Nodes)
		
		
		
		# Find the index of nodes that are fixed
		Index <- NULL
		for (i in 1:length(Fix[[thischain]])){
			var <- names(Fix[[thischain]][i])
				I <- which(All_Nodes_names == var)
				Index <- c(Index,I)
		}
		# Keep only the node that MUST be simulated
		#(Remove the node that are fixed)
		# Nodes_to_simulate <- All_Nodes[-Index]
		
		Nodes_to_simulate <- All_Nodes_names[-Index]
		
		
		
		# Step 2 - Simulate all nodes except those fixed (used as constraint)    -------------------
		# 
		# Simulation for inits 1 ----------------------
		
		# Simulation only of nodes that have to be simulated
		# other nodes are fixed
		set.seed(seed=thischain) 
		mymod$simulate(nodes = Nodes_to_simulate)
		
		
		# head(mymod1)
		# Create inits and save -------------------------------------------
		# 
		
		# Get all stochastic nodes only (only stochastic nodes must really be given a inits)
		# Stochastic nodes that are "data" are also excluded as no inits is required
		# inits_nodes <- mymod$getNodeNames(stochOnly = TRUE, includeData = FALSE, returnType="names")
		# names_inits_nodes <- mymod$getVarNames(includeLogProb = FALSE, nodes=inits_nodes)
		
		inits_nodes <- mymod$getNodeNames(stochOnly = TRUE, includeData = FALSE, returnType="names")
		names_inits_nodes <- mymod$getVarNames(includeLogProb = FALSE, nodes=inits_nodes)
		head(names_inits_nodes)
		
		thisinits <- list()
		
		for (i in 1:length(names_inits_nodes))    {
			thisinits[[i]] <- mymod[[names_inits_nodes[i]]]
		}
		names(thisinits) <- names_inits_nodes
		inits[[thischain]] <- thisinits
		sink()
		sink(type = "message")
		cli_process_done()
		
	}
	# save(inits,file="inits_nimble.RData")
	return(inits)
}







