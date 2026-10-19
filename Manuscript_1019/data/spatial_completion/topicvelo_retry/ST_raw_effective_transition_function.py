def Combined_Topics_Transitions(adata, topics = None, 
                                topic_weights_th_percentile = None,
                                recompute = True,
                                recompute_matrix = True,
                                steady_state_perc = 95, QC_on_topic_genes = False,
                                velocity_type = 'burst',
                                infer_xkey='spliced', infer_ukey='unspliced',
                                embed_xkey = 'Ms', embed_ukey ='Mu',
                                topic_type = 'fastTopics', top_genes_key = 'top_genes', params_key = 'topicVelo_params',
                                transition_matrix_mode = 'count', transition_matrix_name = None,
                                subset_save_prefix = '', save = None):
    '''
    Main function for TopicVelo
    Compute topic-specific transition matrices then integrate them according to topic weights. 

    Construct a global transition matrix from topic transition matrix
    Topic cells are selected. There need to cells overlapping from all the topics. 
   
    Each topic transition matrix is computed using the top topic genes using FastTopics and user-specfied parameters.
    
    The topic transition matrices are combined using cell weights
    (e.g. if a cell is only assigned to topic 1, that cell has weight 1 in topic 1.
    if a cell is assigned to topic 1 with weight 0.6 and topic 2 with weight 0.3, that cell's transition
    matrix will be 2/3*t2 + 1/3*t1 )
    
    The global transition matrix will be row-normalized to 1. 
    
    ***the adata object should contain the top genes from user-specified criteria in adata.uns['top_genes'] 
    
    topics: the topics we want to piece up together 
        default: None, all of the topics will be used
        otherwise, must be a list. A list of one topic is permitted. 
    
    topic_weights_th_percentile: if a list, must be the same dimensions as topics specified (to be implemented. Take a scalar for now)
            
    '''
    #extract the number of cells
    n = adata.n_obs
    #extract the number of topics:
    K = len(adata.uns[top_genes_key])
    #get the top genes
    top_genes = adata.uns[top_genes_key]
    if QC_on_topic_genes:
        try:
            reasonable_top_genes = adata.uns['reasonable_top_genes']
        except:
            print('Need to perform quality control on topic genes')
            return 
    
    #use all topics if topic is none
    use_all_topics = False
    if topics is None:
        use_all_topics = True
        topics = list(range(K))
        
    #add params to adata
    topicVelo_params = {}
    topicVelo_params['topics'] = topics
    topicVelo_params['topic_weights_th_percentile'] = topic_weights_th_percentile 
    topicVelo_params['steady_state_perc'] = steady_state_perc
    topicVelo_params['infer_xkey'] = infer_xkey
    topicVelo_params['infer_ukey'] = infer_ukey
    topicVelo_params['embed_xkey'] = embed_xkey
    topicVelo_params['embed_ukey'] = embed_ukey 
    adata.uns[params_key] = topicVelo_params
    
    #for computing the confidence (coherence within neighborhood of velocity
    topics_cells_velocity_confidence = np.zeros((len(topics),n))
    
    #get topic cells
    ttc_indices, other_cells_indices = get_cells_indices(adata, topics, 
                                                         topic_weights_th_percentile = topic_weights_th_percentile,
                                                         topic_type = topic_type)
    if len(other_cells_indices) > 0:
        adata_other_cells = adata[other_cells_indices,:] 
    #get topic steady state cells
    ttc_ss_indices, transitient_cells_indices = get_cells_indices(adata, topics, 
                                                                  topic_weights_th_percentile = steady_state_perc, 
                                                                  topic_type = topic_type)
    
    #compute with scVelo
    scv.tl.velocity(adata, vkey='velocity')
    
    #store the transition matrices from topics
    TMs = []
    #compute transition matrix for each topic
    for x in range(len(topics)):
        k = topics[x]
        #get the top genes for topic k
        ttg_k = top_genes[k]
        ttg_indices = [adata.var.index.get_loc(gene_name) for gene_name in ttg_k]
        #subset the data and recompute neighbor list in this subset
        adata_subset = adata[ttc_indices[x], ttg_indices]
        scv.pp.neighbors(adata_subset)
        #subset to steady-state cells
        adata_subset_ss = adata[ttc_ss_indices[x], ttg_indices]
        topic_cell_vel_confidence =  np.zeros(n)
        scv.tl.velocity(adata_subset, vkey='velocity')
        
        if QC_on_topic_genes:
            reasonable_genes = reasonable_top_genes[k]
        else:
            reasonable_genes = ttg_k
                
        if velocity_type == 'stochastic':
            #compute velocity, velocity_graph and the transition matrix
            scv.tl.velocity_graph(adata_subset, vkey='velocity', gene_subset=reasonable_genes)
            scv.tl.velocity_confidence(adata_subset, vkey='velocity')
            topic_cell_vel_confidence =  np.zeros(n)
            topic_cell_vel_confidence[ttc_indices[x]] = adata_subset.obs['velocity_confidence']
            topics_cells_velocity_confidence[x] = topic_cell_vel_confidence
            TMs.append(scv.utils.get_transition_matrix(adata_subset, vkey='velocity'))
            
        elif velocity_type == 'burst':
            #burst topic genes and cells. Do not recompute if no need
            save_infer = subset_save_prefix+'T'+str(k)+'_'+infer_xkey+'_'+embed_xkey
            save_path = save_infer + '.npz'
            #recompute if the file does not exist or forced to recompute
            if not exists(save_path) or recompute:
                Burst_Inference(adata_subset_ss, savestring = save_path, report_freq = 50,
                            xkey = infer_xkey, ukey = infer_ukey,
                            vkey = 'burst_velocity')
            inferredParams = np.load(save_path)
            #add the burst_velocity_gamma values
            adata_subset.var['burst_velocity_gamma'] = inferredParams['Optimzal Parameters'][:, 2]
            adata_subset.var['burst_velocity_KLdiv'] = inferredParams['KLdiv']
            if QC_on_topic_genes:
                reasonable_genes = reasonable_top_genes[k]
            else:
                reasonable_genes = ttg_k
            velocity_graph(adata_subset, 
                        transition_matrix_mode=transition_matrix_mode, 
                        xkey = embed_xkey, ukey = embed_ukey,
                        gene_subset = reasonable_genes)  
            
            #add topic cell velocity coherence
            scv.tl.velocity_confidence(adata_subset, vkey='burst_velocity')
            topic_cell_burst_vel_confidence =  np.zeros(n)
            topic_cell_burst_vel_confidence[ttc_indices[x]] = adata_subset.obs['burst_velocity_confidence']
            topics_cells_velocity_confidence[x] = topic_cell_burst_vel_confidence
                
        
            #add topic gene kl divergence (from simulated with MLE to experimental)
            topic_gene_kl_str = 'Topic'+str(k)+'_KLdiv'
            topic_gene_kl = np.zeros(adata.n_vars)
            topic_gene_kl[ttg_indices] = inferredParams['KLdiv']
            adata.var[topic_gene_kl_str] = topic_gene_kl
            #add topic-specific transition 
            TMs.append(scv.utils.get_transition_matrix(adata_subset, vkey='burst_velocity'))

    TM_save_path = subset_save_prefix + 'TransitionMatrix.npz'
    
    #Use topic-specific matrices and topic weights to construct integrated transition matrix
    if not exists(TM_save_path) or recompute_matrix:
        #compute the cell weights for transition matrix
        cells_weights_for_tm = np.zeros((n,len(topics)))
        #iterate through number of topics consider
        for x in range(len(topics)):
            k = topics[x]
            cells_k_indices =  ttc_indices[x]
            cells_weights_k = adata.obs[topic_type+'_'+str(k)]
            #iterate the cells within the topics
            for i in range(len(cells_k_indices)):
                #get the index (in the global adata) of the topic cell
                cell_ki_ind = cells_k_indices[i]
                cells_weights_for_tm[cell_ki_ind, x] = cells_weights_k[cell_ki_ind]
        #row normalize to 1
        for i in range(n):
            weight_sum = np.sum(cells_weights_for_tm[i])
            if weight_sum > 0:
                cells_weights_for_tm[i] = cells_weights_for_tm[i]/weight_sum

        combined_TM = csr_matrix((n,n), dtype=np.float32)
        for x in range(len(topics)):
            cells_k_indices = np.asarray(ttc_indices[x], dtype=int)
            local = TMs[x].tocoo()
            rows = cells_k_indices[local.row]
            cols = cells_k_indices[local.col]
            values = local.data * cells_weights_for_tm[rows, x]
            contribution = csr_matrix((values,(rows,cols)),shape=(n,n))
            combined_TM = (combined_TM + contribution).astype(np.float32)
        sums = np.asarray(combined_TM.sum(axis=1)).ravel()
        empty = np.flatnonzero(sums == 0)
        adata.uns['transition_empty_rows_before_native_fallback'] = empty
        inv = np.divide(1.0,sums,out=np.zeros_like(sums),where=sums!=0)
        combined_TM = combined_TM.multiply(inv[:,None]).tocsr()
        combined_TM.data = np.round(combined_TM.data,6)
        sums = np.asarray(combined_TM.sum(axis=1)).ravel()
        inv = np.divide(1.0,sums,out=np.zeros_like(sums),where=sums!=0)
        combined_TM = combined_TM.multiply(inv[:,None]).tocsr()
        if len(empty):
            conn=adata.obsp['connectivities'].tocsr()
            fallback_rows=[];fallback_cols=[];fallback_values=[]
            for i in empty:
                neighbors=conn.indices[conn.indptr[i]:conn.indptr[i+1]]
                if not len(neighbors): raise ValueError('No fallback neighbors')
                fallback_rows.extend([i]*len(neighbors));fallback_cols.extend(neighbors)
                fallback_values.extend([1/len(neighbors)]*len(neighbors))
            combined_TM += csr_matrix((np.asarray(fallback_values,dtype=np.float32),(fallback_rows,fallback_cols)),shape=(n,n))
        #save transition matrix 
        combined_TM = csr_matrix(combined_TM)
        save_npz(TM_save_path, combined_TM)
        #compute aggregate velocity confidence with cellweights 
        if velocity_type == 'burst':
            adata.obs['topicVelo_velocity_confidence_by_topicWeights'] = np.diag(np.matmul(cells_weights_for_tm, 
                                                                                           topics_cells_velocity_confidence))
            adata.obs['topicVelo_velocity_confidence'] = topics_cells_velocity_confidence.max(axis=0)
            adata.obsm['topicVelo_velocity_confidence_by_topics'] = topics_cells_velocity_confidence.T
        elif velocity_type == 'stochastic':
            adata.obs['scVelo+TM_velocity_confidence_by_topicWeights'] = np.diag(np.matmul(cells_weights_for_tm, 
                                                                                           topics_cells_velocity_confidence))
            adata.obs['scVelo+TM_velocity_confidence'] = topics_cells_velocity_confidence.max(axis=0)
            adata.obsm['scVelo+TM_velocity_confidence_by_topics'] = topics_cells_velocity_confidence.T
    else:
        combined_TM = load_npz(TM_save_path)
        combined_TM = csr_matrix(combined_TM)
    if transition_matrix_name:
        adata.obsp[transition_matrix_name+'_T'] = combined_TM
    elif velocity_type == 'burst':
        adata.obsp['topicVelo_T'] = combined_TM
    elif velocity_type == 'stochastic':
        adata.obsp['topic_modeling+scvelo_stochastic_T'] = combined_TM

    return combined_TM, ttc_indices
