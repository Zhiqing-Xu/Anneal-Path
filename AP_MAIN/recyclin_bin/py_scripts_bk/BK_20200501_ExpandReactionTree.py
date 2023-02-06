# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
    def expand_reaction_tree(self,fwd_smiles_set,fwd_rxn_set,fwd_smiles_set_level_list,trfm_odict,ALL_RXN_SET):
        # Expand the following objects at each level AND UPDATE THE fwd_probabilities_dict:
        # fwd_smiles_set            : intermediate compounds that are ACCEPTED (FOUND and SELECTED), used as reactants at next level
        # fwd_rxn_set             : expands at each level, all reactions FOUND are stored in it
        # processed_rxn_list      : used in reactor.apply_enzymes(), all reactants set are stored in it.
        # fwd_smiles_set_level_list     : all compounds ACCEPTED grouped in different levels, 
        # Could have used the global variables. Those inputs are just for tests.
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        def fwd_SA_probability(hash_a,rxn_list,fwd_smiles_set_level_list,target_cmpd): # Forward Simulated Annealing Probability
            # Return a simulated probability for hash_a, which is a compound found on current level AND UPDATE THE GLOBAL fwd_probabilities_dict
            # 1. Use the following initiated objects to get one similarity score for the reactants (can be many reactants) of one FOUND compound.
            reactants_list=[] # For each reaction contains hash_a as a product, get the reactants.
            simindex_list1=[] # For each reactant in one reaction, get the similarity score
            simindex_list2=[] # For each reaction, get the score of the reactant with the maximum similarity score
            final_reactant_score=0 # Pick the minimum score among all the scores in simindex_list2
            product_score=0
            for i in range(len(rxn_list)):
                if hash_a in rxn_list[i][0]:
                    for j in range(len(rxn_list[i][1])):
                        if rxn_list[i][1][j] in fwd_smiles_set_level_list[subs_level]:
                            reactants_list.append(rxn_list[i][1][j])
                    for one_reactant in reactants_list:
                        simindex_list1.append(self.similarity_score(target_cmpd,one_reactant,fp_type))
                    if simindex_list1!=[]:
                        simindex_list2.append(max(simindex_list1))
                    simindex_list1=[]
                    reactants_list=[]
            if simindex_list2!=[]:
                # ! ! ! ! ! ! ! ! : normally should not be an empty list
                final_reactant_score=max(simindex_list2)
                product_score=self.similarity_score(target_cmpd,hash_a,fp_type)
            else:
                "not normal!!"

            # 2. Calculate probability based on the two scores generated:
            if simindex_list2==[]:
                probability=0
            if product_score>=final_reactant_score:
                probability=1
            else:
                probability=math.exp((product_score-final_reactant_score)/T0_global/(XI_global**subs_level))
            # 3. Update the global fwd_probabilities_dict:
            if fwd_probabilities_dict.has_key(hash_a):
                if probability>fwd_probabilities_dict[hash_a]:
                    fwd_probabilities_dict[hash_a]=probability
            else:
                fwd_probabilities_dict[hash_a]=probability
            return round(probability,3)
	    # ==================================================================================== #
        def fwd_probability_two_cmpds(target_cmpd,reactant,prod):
            # 1. Get similarity score for the reactant and product.
            reactant_score=self.similarity_score(target_cmpd,reactant,fp_type_global)
            prod_score=self.similarity_score(target_cmpd,prod,fp_type_global)
            # 2. Calculate probability based on the two scores generated:
            if prod_score>=reactant_score:
                probability=1
            else:
                probability=math.exp((prod_score-reactant_score)/T0_global/(XI_global**subs_level))
            #print probability, 
            return probability
	    # ==================================================================================== #
        def probability_select(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement
            print "e^(-delta(sim_score)/(T0^subs_level))", "subs_level=", subs_level
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            # 1. option (1) Call fwd_SA_probability to get probabilities for all compounds AND UPDATE THE fwd_probabilities_dict
            '''
            probability_list=[] # [(hash, probability), ( , ), ( , ), ...]
            for cmpd in cmpds_list:
                probability_list.append((cmpd,fwd_SA_probability(cmpd,rxn_list,fwd_smiles_set_level_list,target_cmpd)))
            probability_dict={}
            for (hash,probability) in probability_list:
                probability_dict[hash]=probability
                '''
            # 1. option (2) Update the probability_dict and the global fwd_probabilities_dict directly:
            probability_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                probability=fwd_probability_two_cmpds(target_cmpd,reactant,prod)
                                if (prod in probability_dict.keys() and probability > probability_dict[prod]) or (prod not in probability_dict.keys()):
                                    probability_dict[prod]=probability
                                if (prod in fwd_probabilities_dict.keys() and probability > fwd_probabilities_dict[prod]) or (prod not in fwd_probabilities_dict.keys()):
                                    fwd_probabilities_dict[prod]=probability


            # 2. Rank the probability to accept compounds
            sorted_probability_list = sorted(probability_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print "sorted_probability_list: "
            #print sorted_probability_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_probability_list) else len(sorted_probability_list)
            for j in range(num_cmpds):
                if sorted_probability_list[j][1] != 0:
                    top_similarity_improvement.append(sorted_probability_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement,probability_dict
	    # ==================================================================================== #
        def similarity_improvement_select(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement
            print "e^(-delta(sim_score)/(T0^subs_level))", "subs_level=", subs_level
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            # 1. option (1) Call fwd_SA_probability to get probabilities for all compounds AND UPDATE THE fwd_probabilities_dict
            '''
            probability_list=[] # [(hash, probability), ( , ), ( , ), ...]
            for cmpd in cmpds_list:
                probability_list.append((cmpd,fwd_SA_probability(cmpd,rxn_list,fwd_smiles_set_level_list,target_cmpd)))
            probability_dict={}
            for (hash,probability) in probability_list:
                probability_dict[hash]=probability
                '''
            # 1. option (2) Update the probability_dict and the global fwd_probabilities_dict directly:
            probability_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                probability=fwd_probability_two_cmpds(target_cmpd,reactant,prod)
                                if (prod in probability_dict.keys() and probability > probability_dict[prod]) or (prod not in probability_dict.keys()):
                                    probability_dict[prod]=probability
                                if (prod in fwd_probabilities_dict.keys() and probability > fwd_probabilities_dict[prod]) or (prod not in fwd_probabilities_dict.keys()):
                                    fwd_probabilities_dict[prod]=probability


            # 2. Rank the probability to accept compounds
            sorted_probability_list = sorted(probability_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print "sorted_probability_list: "
            #print sorted_probability_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_probability_list) else len(sorted_probability_list)
            print sorted_probability_list
            for j in range(num_cmpds):
                if sorted_probability_list[j][1] == 1:
                    top_similarity_improvement.append(sorted_probability_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement,probability_dict
	    # ==================================================================================== #
        def similarity_improvement_select_7(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement difference (half number)
            # Used for [-7,-7,-7] alone.
            similarity_improvement_diff_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                similarity_improvement_diff=self.similarity_score(target_cmpd,prod,fp_type_global)-self.similarity_score(target_cmpd,reactant,fp_type_global)

                                if (prod in similarity_improvement_diff_dict.keys() and similarity_improvement_diff > similarity_improvement_diff_dict[prod]) \
                                    or (prod not in similarity_improvement_diff_dict.keys()):
                                    similarity_improvement_diff_dict[prod]=similarity_improvement_diff

            sorted_similarity_improvement_diff_list = sorted(similarity_improvement_diff_dict.items(), key=operator.itemgetter(1), reverse=True)
            print "sorted_similarity_improvement_diff_list: "
            print sorted_similarity_improvement_diff_list
            top_similarity_improvement_diff=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_diff_list) else len(sorted_similarity_improvement_diff_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_diff_list[j][1] != -1:
                    top_similarity_improvement_diff.append(sorted_similarity_improvement_diff_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement_diff
	    # ==================================================================================== #
        def similarity_improvement_select_8(fwd_smiles_set_level_list,cmpds_list,rxn_list,target_cmpd,num_cmpds): # Top ranked similarity improvement ratio (half number)
            # Used for [-8,-8,-8] alone.
            similarity_improvement_ratio_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        for reactant in rxn[1]:
                            if reactant in fwd_smiles_set_level_list[subs_level]:
                                similarity_improvement_ratio=self.similarity_score(target_cmpd,prod,fp_type_global)/(self.similarity_score(target_cmpd,reactant,fp_type_global)+0.00001)

                                if (prod in similarity_improvement_ratio_dict.keys() and similarity_improvement_ratio > similarity_improvement_ratio_dict[prod]) \
                                    or (prod not in similarity_improvement_ratio_dict.keys()):
                                    similarity_improvement_ratio_dict[prod]=similarity_improvement_ratio

            sorted_similarity_improvement_ratio_list = sorted(similarity_improvement_ratio_dict.items(), key=operator.itemgetter(1), reverse=True)
            print "sorted_similarity_improvement_ratio_list: "
            print sorted_similarity_improvement_ratio_list
            top_similarity_improvement_ratio=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_ratio_list) else len(sorted_similarity_improvement_ratio_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_ratio_list[j][1] != 0:
                    top_similarity_improvement_ratio.append(sorted_similarity_improvement_ratio_list[j][0])
                else:
                    break
            # 3. Output list of compounds hashes
            return top_similarity_improvement_ratio
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 0-0 Generate new compounds
        new_rxn_list=self.reactor.apply_enzymes(fwd_smiles_set,fwd_smiles_set_level_list,subs_level,max_C_num,max_O_num)
        #####----------Step 0-1 Update trfm_rxn_dict & trfm_rule_dict.  Modify new_rxn_list. (new version)
        
        modified_new_rxn_set=set([])
        fwd_current_smiles_set=set([])
        fwd_current_accepted_smiles_set=set([])
        for one_rxn in new_rxn_list:
            one_trfm=(one_rxn[0],one_rxn[1])
            if trfm_odict[one_trfm]==0:
                trfm_odict[one_trfm]=set([one_rxn[2],])
                trfm_id=len(trfm_odict)-1
            else:
                trfm_odict[one_trfm].add(one_rxn[2])
                trfm_id=trfm_odict.keys().index(one_trfm)
            one_modified_rxn=(one_rxn[0],one_rxn[1],trfm_id)

            modified_new_rxn_set.add(one_modified_rxn)
            fwd_rxn_set.add(one_modified_rxn)
            fwd_current_smiles_set=fwd_current_smiles_set.union(set(one_modified_rxn[0]))
        new_rxn_list=list(modified_new_rxn_set)
        

        '''
        modified_new_rxn_set=set([])
        fwd_current_smiles_set=set([])
        fwd_current_accepted_smiles_set=set([])
        ALL_RXN_SET=ALL_RXN_SET.union(set(new_rxn_list))
        for one_rxn in new_rxn_list:
            one_modified_rxn=(one_rxn[0],one_rxn[1],"temp_ez_str")
            modified_new_rxn_set.add(one_modified_rxn)
            fwd_rxn_set.add(one_modified_rxn)
            fwd_current_smiles_set=fwd_current_smiles_set.union(set(one_modified_rxn[0]))
        new_rxn_list=list(modified_new_rxn_set)
        '''

        #####----------Step 0-1 (doesn't do anything)
        # use this alternative needs to uncomment step 0-2
        #new_rxn_list=list(set(new_rxn_list))
        #print "len(modified_new_rxn_list): ", len(new_rxn_list)
        print "number of reactions: ", len(new_rxn_list)
        #####----------Step 0-2 Update fwd_rxn_set.  Initialize fwd_current_smiles_set & fwd_current_accepted_smiles_set 
        '''
        fwd_current_smiles_set=set([])
        fwd_current_accepted_smiles_set=set([])
        for one_new_rxn in new_rxn_list:
            fwd_rxn_set.add(one_new_rxn)
            fwd_current_smiles_set=fwd_current_smiles_set.union(set(one_new_rxn[0]))
            '''
        #print "number of compounds found on searching this level (subs side): " + str(len(fwd_current_smiles_set))
        #####----------Step 0-3 Select a number of compounds based on SI and SA (assuming only one target compound here)
        # -1. If last subs level, no need for screening or selecting compounds for next level 
        if subs_level==((max_levels_global+1)/2) -1:
            print "last level, skip compound selection"
            for one_smiles in fwd_current_smiles_set:
                if one_smiles not in fwd_smiles_set:
                    fwd_current_accepted_smiles_set.add(one_smiles)
                    fwd_smiles_set.add(one_smiles)
        else:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
        # 0. Pre-screening: Remove large compounds (have more C or over 2 more O atoms than the target compound)
            screened_smiles_list=[]
            for one_smiles in fwd_current_smiles_set:
                if one_smiles not in bkgd_cmpd_list and one_smiles not in bad_ss_dict.keys():
                    if one_smiles.count("C")<=max_C_num and \
                        one_smiles.count("O")<=max_O_num:
                        screened_smiles_list.append(one_smiles)
            #print "number of compounds screened on searching this level (subs side): " + str(len(screened_smiles_list))
            print "number of compounds: " + str(len(screened_smiles_list))
            if subs_level==0:
                for one_smiles in fwd_current_smiles_set:
                    if one_smiles not in fwd_smiles_set:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)

            else:
                # 0. Determine the number of compounds to accept based on the P0_global value
                if P0_global==[0,0,0]: # Select half number of compounds with higher similarity scores
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[0,0,SISAaccept]
                elif P0_global==[1,1,1]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[2,2,2]: # Select all compounds
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[3,3,3]: # Select based on the probability alone
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[0,SISAaccept,0]
                elif P0_global==[4,4,4]: # Select half number of compounds from three groups (probability group adjusted)
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[-1,-1,-1]: # Select all compounds that improve the similarity (max{products}>min{reactants}, USELESS!)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-6,-6,-6]: # Select all compounds that improve the similarity (max{products}>max{reactants})
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-7,-7,-7]: # Select half number of compounds with higher similarity differences
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-8,-8,-8]: # Select half number of compounds with higher similarity improvement ratios #(***)
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),0,0]
                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[-10,-10,-10]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                else: # Select certain numbers of compounds from the three groups (User-specified)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[min(int(SISAaccept*bin_adj_global[0]),P0_global[0]),\
                                            min(int(SISAaccept*bin_adj_global[1]),P0_global[1]),\
                                            min(int(SISAaccept*bin_adj_global[2]),P0_global[2])]

                
                if P0_global==[0,0,0]: # Select half number of compounds with higher similarity scores
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[1,1,1]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(fwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[2,2,2]: # Select all compounds
                    for one_smiles in screened_smiles_list:
                        if one_smiles not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                    print "Accepted All:", len(screened_smiles_list)

                elif P0_global==[3,3,3]: # (Need to modify this part to make this work)
                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(fwd_probabilities_dict.keys())
                    for one_smiles in random_order_list:
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)

                elif P0_global==[4,4,4]:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)
                    
                    # Probability Select (select according to the probability dict)
                    # The probability_dict used here contains the probabilities of accepting compounds generated by this step.
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(probability_dict.keys())
                    probability_adj_coef=1/sum(probability_dict.values())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if probability_dict[one_smiles]*probability_adj_coef>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[-1,-1,-1]:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=similarity_improvement_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)
                
                elif P0_global==[-6,-6,-6]: # Not edited yet!
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=similarity_improvement_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                elif P0_global==[-7,-7,-7]:
                    # Similarity Improvement ( # Select half number of compounds with higher similarity improvements #(***) )
                    top_similarity_improvement=similarity_improvement_select_7(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                elif P0_global==[-8,-8,-8]:
                    # Similarity Improvement ( # Select half number of compounds with higher similarity improvements #(***) )
                    top_similarity_improvement=similarity_improvement_select_8(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=similarity_improvement_select_7(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    Probabilities_select=[]
                    random_order_list=randomList(fwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[-10,-10,-10]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=similarity_improvement_select_8(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    Probabilities_select=[]
                    random_order_list=randomList(fwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                else:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement,probability_dict)=probability_select(fwd_smiles_set_level_list,screened_smiles_list,new_rxn_list,target_compound,num_screening_global[0])
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        fwd_current_accepted_smiles_set.add(one_smiles)
                        fwd_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(fwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if fwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in fwd_smiles_set:
                            count_x=count_x-1
                            fwd_current_accepted_smiles_set.add(one_smiles)
                            fwd_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,target_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in fwd_smiles_set:
                            fwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            fwd_smiles_set.add(smiles_tuple[0])
                            Similarities_select.append(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        fwd_smiles_set_level_list[subs_level+1]=fwd_current_accepted_smiles_set

        return

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #   
    def bwd_expand_reaction_tree(self,bwd_pathways_list,cmplt_pathways_list,bwd_rxn_list,bwd_level_dict,trfm_odict,ALL_RXN_SET):
        # Expand (perform backward substitution to) the following objects at each level AND UPDATE THE bwd_probabilities_dict:
        # bwd_pathways_list  : shown in the long comments below
        # cmplt_pathways_lis : ???
        # bwd_rxn_list       : expands at each level, all reactions FOUND are stored in it
        # bwd_level_dict     : ???
        # bwd_pathways_list:
            #                [
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme') ],
            #                        [ ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( )  ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #
            #                    [
            #                        [ (-1, ('target_cmpd'), ('reactant_bwd_lv_1_cmpd', ...), 'enzyme'), 
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-2, ('product_bwd_lv_1_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'), ( ), ... ( ), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                                (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-2, ('target_cmpd'), ('reactant_bwd_lv_2_cmpd', ...), 'enzyme'),
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), 
            #                            (-3, ('product_bwd_lv_2_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme'), ( ), ... ( ) ], [ ], ... [ ],
            #                        [ (-3, ('target_cmpd'), ('reactant_bwd_lv_3_cmpd', ...), 'enzyme') ], [ ], ... [ ]
            #                    ],
            #                    ......
            #                ]
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        def bwd_SA_probability(hash_a,rxn_list,starting_cmpds,bwd_level_dict):
            # Return a simulated probability for hash_a, which is a compound found on current level AND UPDATE THE GLOBAL bwd_probabilities_dict
            # 1. Use the following initiated objects to get one similarity score for the reactants (can be many reactants) of one FOUND compound.
            reactants_list=[] # For each reaction contains hash_a as a product, get the reactants.
            simindex_list1=[] # For all reactants in reactants_list, get the similarity score, select the minimum one
            simindex_list2=[] # Get similarity score between hash_a and all starting compounds
            final_reactant_score=0
            product_score=0
            for i in range(len(rxn_list)):
                if hash_a in rxn_list[i][0]:
                    if rxn_list[i][1][0] not in reactants_list: # There shall be only one reactant (for all reactions in rxn_list).
                        reactants_list.append(rxn_list[i][1][0])
            for one_reactant in reactants_list: # Reactants here actually means the products of a reaction in a pathway
                for one_starting_cmpd in starting_cmpds:
                    simindex_list1.append(self.similarity_score(one_starting_cmpd,one_reactant,fp_type_global))
            final_reactant_score=max(simindex_list1)
            for one_starting_cmpd in starting_cmpds:
                simindex_list2.append(self.similarity_score(one_starting_cmpd,hash_a,fp_type_global))
            product_score=max(simindex_list2)
            # 2. Calculate probability based on the two scores generated:
            if product_score>=final_reactant_score:
                probability=1
            else:
                probability=math.exp((product_score-final_reactant_score)/T0_global/(XI_global**prod_level))
            # 3. Update the global bwd_probabilities_dict:
            if bwd_probabilities_dict.has_key(hash_a):
                if probability>bwd_probabilities_dict[hash_a]:
                    bwd_probabilities_dict[hash_a]=probability
            else:
                bwd_probabilities_dict[hash_a]=probability
                #bwd_level_dict[hash_a]=prod_level+1
            return round(probability,3)
        # ==================================================================================== #
        def bwd_probability_two_cmpds(starting_cmpds,reactant,prod):
            # 1. Use the following initiated objects to get similarity score for the reactant and product.
            simindex_list1=[] # All scores for the reactant (reactant tb predicted, 'products' in the pathways)
            simindex_list2=[] # All scores for the prod (prod tb expanded, 'reactants' in the pathways)
            for one_starting_cmpd in starting_cmpds:
                simindex_list1.append(self.similarity_score(one_starting_cmpd,reactant,fp_type_global))
                simindex_list2.append(self.similarity_score(one_starting_cmpd,prod,fp_type_global))
            reactant_score=max(simindex_list1)
            prod_score=max(simindex_list2)
            # 2. Calculate probability based on the two scores generated:
            if prod_score>=reactant_score:
                probability=1
            else:
                probability=math.exp((prod_score-reactant_score)/T0_global/(XI_global**prod_level))
            return probability
        # ==================================================================================== #
        def bwd_probability_select(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            print "e^(-delta(sim_score)/(T0^prod_level))", "prod_level=", prod_level
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            probability_list=[] # [(hash, probability), ( , ), ( , ), ...]
            probability_dict={}
            # 1. option (1) Call bwd_SA_probability to get probabilities for all compounds AND UPDATE THE fwd_probabilities_dict:
            '''
            for cmpd in cmpds_list:
                probability_list.append((cmpd,bwd_SA_probability(cmpd,rxn_list,starting_cmpds,bwd_level_dict)))
            for (hash,probability) in probability_list:
                probability_dict[hash]=probability
                '''

            # 1. option (2) Update the probability_dict and the global bwd_probabilities_dict directly:

            probability_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        probability=bwd_probability_two_cmpds(starting_cmpds,rxn[1][0],prod)
                        if (prod in probability_dict.keys() and probability > probability_dict[prod]) or (prod not in probability_dict.keys()):
                            probability_dict[prod]=probability
                        if (prod in bwd_probabilities_dict.keys() and probability > bwd_probabilities_dict[prod]) or (prod not in bwd_probabilities_dict.keys()):
                            bwd_probabilities_dict[prod]=probability

            # 2. Rank the probability to select compounds
            sorted_probability_list = sorted(probability_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print "sorted_probability_list: "
            #print sorted_probability_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_probability_list) else len(sorted_probability_list)
            for j in range(num_cmpds):
                if sorted_probability_list[j][1] != 0:
                    top_similarity_improvement.append(sorted_probability_list[j][0])
                else:
                    break
            return top_similarity_improvement, probability_dict
        # ==================================================================================== #
        def bwd_similarity_improvement_select(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            print "e^(-delta(sim_score)/(T0^prod_level))", "prod_level=", prod_level
            # Return a list of compounds according to probability rankings, in other words, similarity improvement
            probability_list=[] # [(hash, probability), ( , ), ( , ), ...]
            probability_dict={}
            # 1. option (1) Call bwd_SA_probability to get probabilities for all compounds AND UPDATE THE fwd_probabilities_dict:
            '''
            for cmpd in cmpds_list:
                probability_list.append((cmpd,bwd_SA_probability(cmpd,rxn_list,starting_cmpds,bwd_level_dict)))
            for (hash,probability) in probability_list:
                probability_dict[hash]=probability
                '''

            # 1. option (2) Update the probability_dict and the global bwd_probabilities_dict directly:

            probability_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        probability=bwd_probability_two_cmpds(starting_cmpds,rxn[1][0],prod)
                        if (prod in probability_dict.keys() and probability > probability_dict[prod]) or (prod not in probability_dict.keys()):
                            probability_dict[prod]=probability
                        if (prod in bwd_probabilities_dict.keys() and probability > bwd_probabilities_dict[prod]) or (prod not in bwd_probabilities_dict.keys()):
                            bwd_probabilities_dict[prod]=probability

            # 2. Rank the probability to select compounds
            sorted_probability_list = sorted(probability_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print "sorted_probability_list: "
            #print sorted_probability_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_probability_list) else len(sorted_probability_list)
            print sorted_probability_list
            for j in range(num_cmpds):
                if sorted_probability_list[j][1] ==1:
                    top_similarity_improvement.append(sorted_probability_list[j][0])
                else:
                    break
            return top_similarity_improvement, probability_dict
        # ==================================================================================== #
        def bwd_similarity_improvement_select_7(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            # Used for [-7,-7,-7] alone.
            similarity_improvement_diff_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        similarity_improvement_diff=max([self.similarity_score(prod,one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])- \
                                                        max([self.similarity_score(rxn[1][0],one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])
                        if (prod in similarity_improvement_diff_dict.keys() and similarity_improvement_diff > similarity_improvement_diff_dict[prod]) or (prod not in similarity_improvement_diff_dict.keys()):
                            similarity_improvement_diff_dict[prod]=similarity_improvement_diff
            sorted_similarity_improvement_diff_list = sorted(similarity_improvement_diff_dict.items(), key=operator.itemgetter(1), reverse=True)
            print "sorted_similarity_improvement_diff_list: "
            print sorted_similarity_improvement_diff_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_diff_list) else len(sorted_similarity_improvement_diff_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_diff_list[j][1] !=-1:
                    top_similarity_improvement.append(sorted_similarity_improvement_diff_list[j][0])
                else:
                    break
            return top_similarity_improvement
        # ==================================================================================== #
        def bwd_similarity_improvement_select_8(cmpds_list,rxn_list,starting_cmpds,num_cmpds,bwd_level_dict):
            # Used for [-8,-8,-8] alone.
            similarity_improvement_ratio_dict={}
            for rxn in rxn_list:
                for prod in rxn[0]:
                    if prod in cmpds_list:
                        similarity_improvement_ratio=max([self.similarity_score(prod,one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds])/ \
                                                        max([self.similarity_score(rxn[1][0],one_starting_cmpd,fp_type_global) for one_starting_cmpd in starting_cmpds]+[0.00001])
                        if (prod in similarity_improvement_ratio_dict.keys() and similarity_improvement_ratio > similarity_improvement_ratio_dict[prod]) or (prod not in similarity_improvement_ratio_dict.keys()):
                            similarity_improvement_ratio_dict[prod]=similarity_improvement_ratio
            sorted_similarity_improvement_ratio_list = sorted(similarity_improvement_ratio_dict.items(), key=operator.itemgetter(1), reverse=True)
            print "sorted_similarity_improvement_ratio_list: "
            print sorted_similarity_improvement_ratio_list
            top_similarity_improvement=[] # Output list.
            num_cmpds= int(num_cmpds) if num_cmpds<len(sorted_similarity_improvement_ratio_list) else len(sorted_similarity_improvement_ratio_list)
            for j in range(num_cmpds):
                if sorted_similarity_improvement_ratio_list[j][1] !=0:
                    top_similarity_improvement.append(sorted_similarity_improvement_ratio_list[j][0])
                else:
                    break
            return top_similarity_improvement
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 0-0 Draw all reactants to be substitued and prepare for backward reaction search
        cmpds_need_expand=[] # A list of compounds input into reactor.bwd_apply_enzymes()
        for i in range(len(bwd_pathways_list[-1])):
            for j in range(len(bwd_pathways_list[-1][i])):
                if bwd_pathways_list[-1][i][j][0] == 0 - prod_level:
                    for k in range(len(bwd_pathways_list[-1][i][j][2])):
                        if bwd_pathways_list[-1][i][j][2][k] not in subs_smiles and bwd_pathways_list[-1][i][j][2][k] not in bkgd_cmpd_list:
                            cmpds_need_expand.append(bwd_pathways_list[-1][i][j][2][k])
        screened_cmpds_need_expand_list=[]
        for one_smiles in cmpds_need_expand:
            if one_smiles not in bkgd_cmpd_list:
                if one_smiles.count("C")<=max_C_num and \
                    one_smiles.count("O")<=max_O_num:
                    screened_cmpds_need_expand_list.append(one_smiles)

        #####----------Step 0-1 Generate new compounds and update reaction list
        # bwd_new_rxn_list : A list of compounds contains outputs of reactor.bwd_apply_enzymes()
        bwd_new_rxn_list=self.reactor.bwd_apply_enzymes_zx(screened_cmpds_need_expand_list,prod_level)

        #####----------Step 0-2 Update trfm_rxn_dict & trfm_rule_dict.  Modify bwd_new_rxn_list. (new version)
        
        modified_bwd_new_rxn_set=set([])
        for one_rxn in bwd_new_rxn_list:
            one_trfm=(one_rxn[0],one_rxn[1])
            if trfm_odict[one_trfm]==0:
                trfm_odict[one_trfm]=set([one_rxn[2],])
                trfm_id=len(trfm_odict)-1
            else:
                trfm_odict[one_trfm].add(one_rxn[2])
                trfm_id=trfm_odict.keys().index(one_trfm)
            modified_bwd_new_rxn_set.add((one_rxn[0],one_rxn[1],trfm_id))
        bwd_new_rxn_list=list(modified_bwd_new_rxn_set)
        
        '''
        modified_bwd_new_rxn_set=set([])
        ALL_RXN_SET=ALL_RXN_SET.union(set(bwd_new_rxn_list))
        for one_rxn in bwd_new_rxn_list:
            modified_bwd_new_rxn_set.add((one_rxn[0],one_rxn[1],"temp_ez_str"))
        bwd_new_rxn_list=list(modified_bwd_new_rxn_set)
        '''
        #####----------
        #bwd_new_rxn_list=list(set(bwd_new_rxn_list))


        #####----------Step 0-2 Update bwd_rxn_set.  Initialize bwd_current_smiles_set & bwd_current_accepted_smiles_set 
        bwd_current_smiles_set=set([])
        bwd_current_accepted_smiles_set=set([])

        for i in range(len(bwd_new_rxn_list)):
            bwd_rxn_list.append(bwd_new_rxn_list[i])
            bwd_current_smiles_set=bwd_current_smiles_set.union(set(bwd_new_rxn_list[i][0]))

        #####----------Step 0-3 Draw all products from predicted reactions and get ready for SI and SA selection
        #print "number of compounds found on searching this level (prod side): " + str(len(bwd_current_smiles_set))
        for one_smiles in bwd_current_smiles_set:
            if one_smiles not in bwd_level_dict.keys():
                bwd_level_dict[one_smiles]=prod_level+1
        #####----------Step 0-4 Select (accept) a number of compounds based on SI and SA
        # -1. If last prod level, no need for screening or selecting compounds for next level 
        if prod_level==((max_levels_global)/2) -1:
            print "last level, skip compound selection"
            bwd_current_accepted_smiles_set=bwd_current_accepted_smiles_set.union(bwd_current_smiles_set)
        else:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
        # 0. Pre-screening: actually don't need to remove large compounds for backward search! So could have remove this step here!
            screened_smiles_list=[]
            for one_smiles in bwd_current_smiles_set:
                if one_smiles not in bkgd_cmpd_list:
                    if one_smiles.count("C")<=max_C_num and \
                        one_smiles.count("O")<=max_O_num:
                        screened_smiles_list.append(one_smiles)

            #print "number of compounds screened on searching this level (subs side): " + str(len(screened_smiles_list))
            print "number of compounds: " + str(len(screened_smiles_list))
            if prod_level==0:
                for one_smiles in bwd_current_smiles_set:
                    bwd_current_accepted_smiles_set.add(one_smiles)
            else:
                # 0. Determine the number of compounds to accept based on the P0_global value
                if P0_global==[0,0,0]: # Select half number of compounds with higher similarity scores
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[0,0,SISAaccept]
                elif P0_global==[1,1,1]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[2,2,2]: # Select all compounds
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[3,3,3]: # Select based on the probability alone
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[0,SISAaccept,0]
                elif P0_global==[4,4,4]: # Select half number of compounds from three groups (probability group adjusted)
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[-1,-1,-1]: # Select all compounds that improve the similarity (max{products}>min{reactants}, USELESS!)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-6,-6,-6]: # Select all compounds that improve the similarity (max{products}>max{reactants})
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-7,-7,-7]: # Select half number of compounds with higher similarity differences
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[SISAaccept,0,0]
                elif P0_global==[-8,-8,-8]: # Select half number of compounds with higher similarity improvement ratios #(***)
                    SISAaccept=(len(screened_smiles_list)/2)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),0,0]
                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                elif P0_global==[-10,-10,-10]: # Select half number of compounds from three groups
                    SISAaccept=(len(screened_smiles_list)/6)+1
                    num_screening_global=[int(SISAaccept*bin_adj_global[0]),int(SISAaccept*bin_adj_global[1]),int(SISAaccept*bin_adj_global[2])]
                else: # Select certain numbers of compounds from the three groups (User-specified)
                    SISAaccept=(len(screened_smiles_list))
                    num_screening_global=[min(int(SISAaccept*bin_adj_global[0]),P0_global[0]),\
                                            min(int(SISAaccept*bin_adj_global[1]),P0_global[1]),\
                                            min(int(SISAaccept*bin_adj_global[2]),P0_global[2])]

                
                if P0_global==[0,0,0]: # Select half number of compounds with higher similarity scores
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[1,1,1]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(bwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[2,2,2]: # Select all compounds
                    bwd_current_accepted_smiles_set=screened_smiles_list
                    print "Accepted All:", 
                    print len(screened_smiles_list)

                elif P0_global==[3,3,3]:
                    # Probability Select (select according to the probability dict)
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(bwd_probabilities_dict.keys())
                    for one_smiles in random_order_list:
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)

                elif P0_global==[4,4,4]:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)
                    
                    # Probability Select (select according to the probability dict)
                    # The probability_dict used here contains the probabilities of accepting compounds generated by this step.
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(probability_dict.keys())
                    probability_adj_coef=1/sum(probability_dict.values())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if probability_dict[one_smiles]*probability_adj_coef>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[-1,-1,-1]:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement, probability_dict)=bwd_similarity_improvement_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)

                elif P0_global==[-7,-7,-7]:
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_7(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                elif P0_global==[-8,-8,-8]:
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_8(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                elif P0_global==[-9,-9,-9]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_7(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    print "Probabilities: "
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    Probabilities_select=[]
                    random_order_list=randomList(bwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                elif P0_global==[-10,-10,-10]: # Select half number of compounds from three groups
                    # Similarity Improvement (select according to the probability dict)
                    top_similarity_improvement=bwd_similarity_improvement_select_8(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                    print len(top_similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    print "Probabilities: "
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    Probabilities_select=[]
                    random_order_list=randomList(bwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)

                else:
                    # Similarity Improvement (select according to the probability dict)
                    (top_similarity_improvement, probability_dict)=bwd_probability_select(screened_smiles_list,bwd_new_rxn_list,subs_smiles,num_screening_global[0],bwd_level_dict)
                    print "Similarity improvement: "
                    Similarity_improvement=[]
                    for one_smiles in top_similarity_improvement:
                        bwd_current_accepted_smiles_set.add(one_smiles)
                        Similarity_improvement.append(one_smiles)
                    print len(Similarity_improvement)

                    # Probability Select (select according to the probability dict)
                    # The bwd_probabilities_dict used here contains the probabilities of accepting compounds (ALL generated).
                    print "Probabilities: "
                    Probabilities_select=[]
                    random_order_list=randomList(bwd_probabilities_dict.keys())
                    count_x=num_screening_global[1]
                    for one_smiles in random_order_list*10:
                        if count_x==0:
                            break
                        if bwd_probabilities_dict[one_smiles]>=random.random() and one_smiles not in bwd_current_accepted_smiles_set:
                            count_x=count_x-1
                            bwd_current_accepted_smiles_set.add(one_smiles)
                            Probabilities_select.append(one_smiles)
                    print len(Probabilities_select)
                    
                    # Similarity Select (based on similarity scores)
                    print "Similarities: "
                    Similarities_select=[]
                    taofactors_dict=self.similarity_dict(screened_smiles_list,subs_smiles,fp_type_global,num_screening_global[2])
                    sorted_taofactors_list = sorted(taofactors_dict.items(), key=operator.itemgetter(1), reverse=True)
                    top_ranked_reactants=[]
                    num_cmpds= num_screening_global[2] if num_screening_global[2]<len(sorted_taofactors_list) else len(sorted_taofactors_list)
                    for smiles_tuple in sorted_taofactors_list:
                        if smiles_tuple[0] not in bwd_current_accepted_smiles_set:
                            bwd_current_accepted_smiles_set.add(smiles_tuple[0])
                            num_cmpds=num_cmpds-1
                            Similarities_select.append(smiles_tuple[0])
                            #print smiles_tuple[0]
                        if num_cmpds<=0:
                            break
                    print len(Similarities_select)
        
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
        #####----------Step 1 Take the newest group of incomplete pathways from bwd_pathways_list (in last [])
        tb_expanded_list=bwd_pathways_list[-1] # get the last list in bwd_pathways_list, that is the list to be expanded upon
        bwd_pathways_list.append([]) # Prepare for adding a new group of incomplete pathways of current level

        #####----------Step 2 Expand 2 pathway lists for backward direction search(, FIRST add reactions with compounds found on THIS LEVEL).
        if bwd_pathways_list[0][0][0][0]==0:
            #####----------Step 2-0 
            # At the first backward search step, initiate the bwd_pathways_list with the first backward level search result
            bwd_pathways_list=[[]] # 
            for n in range(len(bwd_new_rxn_list)):
                for o in range(len(bwd_new_rxn_list[n][0])):
                    if bwd_new_rxn_list[n][0][o] in bwd_current_accepted_smiles_set:
                        one_pathway_list=[]
                        one_pathway_list.append((-1,bwd_new_rxn_list[n][1],bwd_new_rxn_list[n][0],bwd_new_rxn_list[n][2]))
                        bwd_pathways_list[-1].append(one_pathway_list)
                        del one_pathway_list
                        break

        else:
            #####----------Step 2-1 For both full-length and non-full-length pathways, add reactions with compounds found on THIS LEVEL.
            # After the first backward search step, bwd_pathways_list now has been initiated and are ready for backward substitution.
            print"Compounds to be expanded", len(tb_expanded_list)

            for n in range(len(tb_expanded_list)):
                # 0. Obtain all compounds along this one certain pathway, store in reactants_found_already.
                reactants_found_already=set([])
                for o in range(len(tb_expanded_list[n])):
                    for p in range(len(tb_expanded_list[n][o])):
                        if p==1 or p==2:
                            reactants_found_already = reactants_found_already | set(tb_expanded_list[n][o][p])

                # 1. Obtain all compounds to be substituted by predicted reactants(, duplicates and starting compounds removed).
                reactantstobeexpanded=[]
                for o in range(len(tb_expanded_list[n])):
                    if tb_expanded_list[n][o][0] == 0 - prod_level:
                        for p in range(len(tb_expanded_list[n][o][2])):
                            if (tb_expanded_list[n][o][2][p] not in reactantstobeexpanded) and tb_expanded_list[n][o][2][p] not in subs_smiles:
                                reactantstobeexpanded.append(tb_expanded_list[n][o][2][p])

                # 2. Normally reactantstobeexpanded shall not be empty list.
                if len(reactantstobeexpanded)==0:
                    print "completed pathway found in incomplete pathway list"
                    continue

                # 3. For each compounds to be substituted, get all reactions (from this step of search) that can be used to substitute this compound.
                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_new_rxn_list)):
                        for s in range(len(bwd_new_rxn_list[r][0])):
                            if bwd_new_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_new_rxn_list[r][1] and \
                                bwd_new_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_new_rxn_list[r][0][s]] == prod_level+1:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on this level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_new_rxn_list[r][1],bwd_new_rxn_list[r][0],bwd_new_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.
                



                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_new_rxn_list)):
                        for s in range(len(bwd_new_rxn_list[r][0])):
                            if bwd_new_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_new_rxn_list[r][1] and \
                                bwd_new_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_new_rxn_list[r][0][s]] == prod_level+1:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on this level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_new_rxn_list[r][1],bwd_new_rxn_list[r][0],bwd_new_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.



                # 4. Use cart_prod to get all combinations of those reactions that is all possibilities how compounds can be substituted
                next_level_reactions_set=cart_prod(next_level_reactions)
                # 5. Expand one pathway and check if it is completed
                for u in range(len(next_level_reactions_set)):
                    one_pathway_list=deepcopy(tb_expanded_list[n])
                    for v in range(len(next_level_reactions_set[u])):
                        if next_level_reactions_set[u][v] not in one_pathway_list: # Prevent adding two identical reactions into one pathway (, doesnt do anything to avoid bad pathways)
                            one_pathway_list.append(next_level_reactions_set[u][v])
                    if one_pathway_list not in bwd_pathways_list[prod_level] and one_pathway_list not in cmplt_pathways_list:
                        # 5.1 Now check if the pathway is already completed.
                        reactantswillbeexpanded=[] # Predict all compounds to be substituted at next level, if none, pathway is completed.
                        for w in range(len(one_pathway_list)):
                            if one_pathway_list[w][0] == 0 - prod_level - 1:
                                for x in range(len(one_pathway_list[w][2])):
                                    if (one_pathway_list[w][2][x] not in reactantswillbeexpanded) and one_pathway_list[w][2][x] not in subs_smiles:
                                        reactantswillbeexpanded.append(one_pathway_list[w][2][x])
                        # 5.2 If a complete pathway, add to another list not the bwd_pahways_list
                        if reactantswillbeexpanded!=[]:
                            bwd_pathways_list[prod_level].append(one_pathway_list)
                        else:
                            cmplt_pathways_list.append(one_pathway_list)
            

        # ! ! ! ! ! ! ! ! : Step 3 has almost been rewritten, so may have mistakes unfound.
        #####----------Step 3 Expand 2 pathway lists for backward direction search (, NOW add reactions with compounds found on PREVIOUS LEVELS).

        for z in range(prod_level):
            print "prod_level", z
            cmpd_level=z+1 # cmpd_level is the level of the compounds to be substituted.
            #####----------Step 3-0 Compounds FOUND at FIRST step get DROPPED and then ACCEPTED at this current level
            if cmpd_level==1: # First deal with compounds FOUND at FIRST step (, this means another new pathway needs to be initiated).
                for n in range(len(bwd_rxn_list)):
                    for o in range(len(bwd_rxn_list[n][0])):
                        if (bwd_rxn_list[n][0][o] in bwd_current_accepted_smiles_set and \
                            bwd_level_dict[bwd_rxn_list[n][0][o]]==1 and \
                            bwd_rxn_list[n][1][0]==target_compound): # ! ! ! ! ! ! ! ! : This line makes sure only compounds found on FIRST level is substituted.
                            one_pathway_list=[]
                            one_pathway_list.append((-1-prod_level,bwd_rxn_list[n][1],bwd_rxn_list[n][0],bwd_rxn_list[n][2]))
                            bwd_pathways_list[-1].append(one_pathway_list)
                            break
                continue # Use a continue here so that dont need a long else part below.
            #print("done step #1.-2")
            #####----------Step 3-1 Compounds FOUND AFTER 1st step and before current step get DROPPED and then ACCEPTED at this current level
            # Prepare for expanding the bwd_pathways_list:
            tb_expanded_list=[]
            for m in range(len(bwd_pathways_list[-1-prod_level+cmpd_level-2])):
                # bwd_pathways_list contains all incomplete pathways histories at each level, stored in a different list.
                # Need to go to one earlier incomplete pathways list, so that subsitution can be made to add compounds found earlier but just accepted.
                # Need to adjust the earlier incomplete pathways to make substitution.
                adjusted_earlier_pathway=[]
                for n in range(len(bwd_pathways_list[-1-prod_level+cmpd_level-2][m])):
                    # adjust the reaction level index to prepare for substitution (use a,b below):
                    a=list(bwd_pathways_list[-1-prod_level+cmpd_level-2][m][n])
                    a[0]=a[0]-prod_level-1+cmpd_level
                    b=tuple(a)
                    adjusted_earlier_pathway.append(b)
                tb_expanded_list.append(adjusted_earlier_pathway)
            #print("done step #1.-1")

            #print "len(tb_expanded_list)", len(tb_expanded_list)
            for n in range(len(tb_expanded_list)):
                #print n
                # 0. Obtain all compounds along this one certain pathway, store in reactants_found_already.
                reactants_found_already=set([])
                for o in range(len(tb_expanded_list[n])):
                    for p in range(len(tb_expanded_list[n][o])):
                        if p==1 or p==2:
                            reactants_found_already = reactants_found_already | set(tb_expanded_list[n][o][p])

                # 1. Obtain all compounds to be substituted by predicted reactants(, duplicates and starting compounds removed).
                reactantstobeexpanded=[]
                for o in range(len(tb_expanded_list[n])):
                    if tb_expanded_list[n][o][0] == 0 - prod_level:
                        for p in range(len(tb_expanded_list[n][o][2])):
                            if (tb_expanded_list[n][o][2][p] not in reactantstobeexpanded) and tb_expanded_list[n][o][2][p] not in subs_smiles:
                                reactantstobeexpanded.append(tb_expanded_list[n][o][2][p])
                # 2. Normally reactantstobeexpanded shall not be empty list.
                if len(reactantstobeexpanded)==0:
                    print "completed pathway found in incomplete pathway list2: "
                    #print tb_expanded_list[n]
                    continue
                # 3. For each compounds to be substituted, get all reactions (from this step of search) that can be used to substitute this compound.
                next_level_reactions=[]
                for q in range(len(reactantstobeexpanded)):
                    next_level_reactions.append([]) # Each reactants has a corresponding list in next_level_reactions, which contains all reactions that can be used to substitute this reactant.
                    for r in range(len(bwd_rxn_list)):
                        for s in range(len(bwd_rxn_list[r][0])):
                            if bwd_rxn_list[r][0][s] in bwd_current_accepted_smiles_set and \
                                reactantstobeexpanded[q] in bwd_rxn_list[r][1] and \
                                bwd_rxn_list[r][0][s] not in reactants_found_already and \
                                bwd_level_dict[bwd_rxn_list[r][0][s]] == cmpd_level:# ! ! ! ! ! ! ! ! : This line makes sure only compounds found on that earlier level is substituted.
                                next_level_reactions[q].append((0 - prod_level -1, bwd_rxn_list[r][1],bwd_rxn_list[r][0],bwd_rxn_list[r][2]))
                                break # Once a reaction is accepted, the other reactants (, or backwardly products) dont need to be checked.
                # 4. Use cart_prod to get all combinations of those reactions that is all possibilities how compounds can be substituted
                next_level_reactions_set=cart_prod(next_level_reactions)
                # 5. Expand one pathway and check if it is completed
                for u in range(len(next_level_reactions_set)):
                    one_pathway_list=deepcopy(tb_expanded_list[n])
                    for v in range(len(next_level_reactions_set[u])):
                        if next_level_reactions_set[u][v] not in one_pathway_list:
                            one_pathway_list.append(next_level_reactions_set[u][v])
                    if one_pathway_list not in bwd_pathways_list[prod_level] and one_pathway_list not in cmplt_pathways_list:
                        # 5.1 Now check if the pathway is already completed.
                        reactantswillbeexpanded=[] # Predict all compounds to be substituted at next level, if none, pathway is completed.
                        for w in range(len(one_pathway_list)):
                            if one_pathway_list[w][0] == 0 - prod_level - 1:
                                for x in range(len(one_pathway_list[w][2])):
                                    if (one_pathway_list[w][2][x] not in reactantswillbeexpanded) and one_pathway_list[w][2][x] not in subs_smiles:
                                        reactantswillbeexpanded.append(one_pathway_list[w][2][x])
                        # 5.2 If a complete pathway, add to another list not the bwd_pahways_list
                        if reactantswillbeexpanded!=[]:
                            bwd_pathways_list[prod_level].append(one_pathway_list)
                        else:
                            cmplt_pathways_list.append(one_pathway_list)

        return bwd_pathways_list
