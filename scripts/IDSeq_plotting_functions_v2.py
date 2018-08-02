import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns

def per_sample_result_plots(input_matrix, filename, nrow, PROBABILITY_THRESHOLD, graph_axis, label, positive_class):
    ''' Split the predictions by patient and plot independent panels of prediction '''
    

    fig, axarr = plt.subplots(nrow,5, figsize=(20, 4.4*nrow), sharex=True, sharey=True)#, facecolor='white')
    row = 0
    col = 0

    c = ['blue','red']
    c3 = ['white','red']

    for pt in list(set(input_matrix['sampleID'])):

        if row < 21 and col < 5:

            #plt.figure(figsize=[8,8])

            curr_plot_sample = input_matrix[input_matrix['sampleID'] == pt]

            # visualize the test set in RNA rM by DNA rM space;
            axarr[row,col].scatter(curr_plot_sample[graph_axis],
                        curr_plot_sample['score'],
                        edgecolor=[c[int(curr_plot_sample.loc[i]['positive'])] 
                                   if curr_plot_sample[curr_plot_sample['sampleID'] == 
                                                  curr_plot_sample.loc[i]['sampleID']]['positive'].sum() < 2 
                                   or c[int(curr_plot_sample.loc[i]['positive'])] == 'blue' 
                                   or curr_plot_sample.loc[i]['score'] == max(
                                       curr_plot_sample[(curr_plot_sample.sampleID == curr_plot_sample.loc[i]['sampleID']) 
                                                   & (curr_plot_sample.positive == True)]['score']) 
                                   else 'orange' for i in curr_plot_sample.index],
                        s=[max(i*600,3) for i in curr_plot_sample['score']],
                        facecolor=[c[int(curr_plot_sample['positive'].iloc[i])] if curr_plot_sample['score'].iloc[i] > PROBABILITY_THRESHOLD else 'white' for i in range(len(curr_plot_sample['positive'])) ],   
                        alpha=.7) 

            for i in curr_plot_sample.index:
                if curr_plot_sample['positive'].loc[i] == True or curr_plot_sample['score'].loc[i] > PROBABILITY_THRESHOLD:
                    sp = curr_plot_sample.loc[i]['name'].split(' ')
                    sn = curr_plot_sample.loc[i]['sampleID'].split(' ')[0]
                    axarr[row,col].annotate(curr_plot_sample.loc[i]['name'], #sp[0][0] + '. ' + sp[1],
                                     (curr_plot_sample.loc[i][graph_axis],curr_plot_sample.loc[i]['score'] + .05), 
                                     fontsize=7, color='grey',rotation=15)     


            axarr[row,col].set_title(str(pt + "\n" + list(set(curr_plot_sample['known_organisms']))[0]), color=c[list(set(curr_plot_sample[label]))[0]== positive_class])
            #axarr[row,col].set_title(pt, color=c[list(set(curr_plot_sample['LRTIstatus']))[0]=='LRTI+C+M'])
            axarr[row,col].set_xticks(np.arange(0,5.5,1))
            axarr[row,col].set_yticks(np.arange(0,1,.2))

            if col == 4:
                row += 1
                col = 0
            else:
                col += 1

    plt.savefig(filename)


def per_sample_result_plots_sorted(input_matrix, filename, nrow, PROBABILITY_THRESHOLD, graph_axis, label, positive_class, negative_class, unknown_class):
    ''' Split the predictions by patient and plot independent panels of prediction '''
    

    fig, axarr = plt.subplots(nrow,5, figsize=(20, 4.4*nrow), sharex=True, sharey=True)#, facecolor='white')
    row = 0
    col = 0

    c = ['blue','red']
    c3 = ['white','red']

    color_count = {}
    color_count[positive_class] = 'red' 
    color_count[negative_class] = 'blue'
    for i in unknown_class:  
      color_count[i] = 'gold'

    input_matrix = input_matrix.sort_values(by='sampleID')

    for classID in list(set(input_matrix[label])):

      for pt in list(set(input_matrix['sampleID'])):

          if row < 21 and col < 5:

              #plt.figure(figsize=[8,8])

              curr_plot_sample = input_matrix[input_matrix['sampleID'] == pt]

              if list(curr_plot_sample[label])[0] == classID:

                # visualize the test set in RNA rM by DNA rM space;
                axarr[row,col].scatter(curr_plot_sample[graph_axis],
                            curr_plot_sample['score'],
                            edgecolor=[c[int(curr_plot_sample.loc[i]['positive'])] 
                                       if curr_plot_sample[curr_plot_sample['sampleID'] == 
                                                      curr_plot_sample.loc[i]['sampleID']]['positive'].sum() < 2 
                                       or c[int(curr_plot_sample.loc[i]['positive'])] == 'blue' 
                                       or curr_plot_sample.loc[i]['score'] == max(
                                           curr_plot_sample[(curr_plot_sample.sampleID == curr_plot_sample.loc[i]['sampleID']) 
                                                       & (curr_plot_sample.positive == True)]['score']) 
                                       else 'orange' for i in curr_plot_sample.index],
                            s=[max(i*600,3) for i in curr_plot_sample['score']],
                            facecolor=[c[int(curr_plot_sample['positive'].iloc[i])] if curr_plot_sample['score'].iloc[i] > PROBABILITY_THRESHOLD else 'white' for i in range(len(curr_plot_sample['positive'])) ],   
                            alpha=.7) 

                for i in curr_plot_sample.index:
                    if curr_plot_sample['positive'].loc[i] == True or curr_plot_sample['score'].loc[i] > PROBABILITY_THRESHOLD:
                        sp = curr_plot_sample.loc[i]['name'].split(' ')
                        sn = curr_plot_sample.loc[i]['sampleID'].split(' ')[0]
                        axarr[row,col].annotate(curr_plot_sample.loc[i]['name'], #sp[0][0] + '. ' + sp[1],
                                         (curr_plot_sample.loc[i][graph_axis],curr_plot_sample.loc[i]['score'] + .05), 
                                         fontsize=7, color='grey',rotation=15)     


                axarr[row,col].set_title(str(pt + "\n" + list(set(curr_plot_sample['known_organisms']))[0]), color=color_count[list(set(curr_plot_sample[label]))[0]])
                #axarr[row,col].set_title(pt, color=c[list(set(curr_plot_sample['LRTIstatus']))[0]=='LRTI+C+M'])
                axarr[row,col].set_xticks(np.arange(0,5.5,1))
                axarr[row,col].set_yticks(np.arange(0,1,.2))

                if col == 4:
                    row += 1
                    col = 0
                else:
                    col += 1

    plt.savefig(filename)