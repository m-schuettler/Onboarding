import gc  # garbage collector
import glob
import matplotlib.pyplot as plt
import napari
import numpy as np
import os
import pandas as pd
from scipy import ndimage as ndi
import skimage as ski
from skimage import io




def count_nuclei_v3(image, save=False, showplots=True):
## version 3: for segmenting DAPI images
# arguments:    'image' should be the filepath to the desired image
#               'save=False' suppressed saving of diagnostic figure
#               'showplots=False' suppresses showing of diagnostic figure when running the function in Jupyter Notebook
    img = io.imread(image)
    grayimg = io.imread(image, as_gray=True)                                    # create image for diagnostic figure
    grayimg = ski.exposure.rescale_intensity(grayimg)
    preproc = ski.filters.gaussian(img, sigma=2.5, preserve_range=True)         # Gaussian blur reduces heterogeneous intensity in nuclei
    preproc = ski.exposure.rescale_intensity(preproc)                           # improves contrast

    # subtract background
    bg = ski.filters.gaussian(preproc, sigma=200)
    preproc = preproc - bg
    preproc[preproc < 0]  = 0  # make sure there are no negative values
    preproc = preproc / np.max(preproc)  # normalize

    # segmentation
    threshold = ski.filters.threshold_triangle(preproc)                         # segmentation with triangle
    threshold = threshold*1.25                                                  # reduce threshold for more detailed segmentation
    seg = ski.morphology.closing(preproc > threshold, ski.morphology.square(3))
    seg = ski.morphology.opening(seg)
    #seg = ski.morphology.remove_small_objects(seg, min_size=100)               # remove small objects under min_size pixels
    
    # watershed
    distance = ndi.distance_transform_edt(seg)                                  # watershed bw to separate touching nuclei
    coords = ski.feature.peak_local_max(distance, footprint=np.ones((4, 4)), labels=seg, min_distance=20)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    wsh = ski.segmentation.watershed(-distance, markers, mask=seg)
    wsh = ski.morphology.remove_small_objects(wsh, min_size=150)                # remove small objects under min_size pixels
    largeobjs = ski.morphology.remove_small_objects(wsh, min_size=3500)         # remove large objects above min_size pixels
    wsh = wsh ^ largeobjs

    _, num = ski.measure.label(wsh, return_num=True)                            # label segmented regions

    # draw coordinates in plot
    name = os.path.basename(image).removeprefix('out_opt_flow_registered_').removesuffix('.tif')
    plt.title(name + '   |   nuclei found: ' + str(num))

    plt.imshow(ski.color.label2rgb(wsh, grayimg, alpha=0.25, bg_label=0, bg_color=None))

    if save:  # if 'save' is 'True'
        if not os.path.exists('../results/Task-14/'):
            os.makedirs('../results/Task-14/')
        plt.savefig('../results/Task-14/Diagnostic-Figure_' + name + '.png')

    if showplots:
        plt.show()
    else:
        plt.close()

    # return results
    return num




def batch_count_nuclei(images):
## applies 'count_nuclei_v3()' to all files in 'images'
# arguments:    'images' should be a list that contains all desired filenames
    # get nucleus counts and names
    numlist = []
    nameslist = []
    for i in range(0, len(images)):
        num = count_nuclei_v3(images[i], save=True, showplots=False)
        numlist.append(num)
        nameslist.append(os.path.basename(images[i]).removeprefix('out_opt_flow_registered_').removesuffix('.tif'))
    
    # put nucleus counts and names together into dataframe
    ziplist = list(zip(nameslist, numlist))
    df = pd.DataFrame(ziplist, columns = ['files', 'nucleus counts'])
    df.to_csv('../results/Task-14/Nucleus-Counts.csv')
    
    return df




def count_nuclei_v4(image, labelspath, figpath, showplots=True):
## version 4:   contour segmentation
##              differenciate segmentation workflow for Atto dye
# arguments:    'image' is the filepath to the desired image
#               'subfolder' is defined in 'bashdebug.sh' and is the name of the subfolder that results are saved to
#               'showplots=False' prevents the figures from being shown
    img = io.imread(image, dtype=np.float64)
    grayimg = io.imread(image, as_gray=True, dtype=np.float64)                          # create image for diagnostic figure
    grayimg = ski.exposure.rescale_intensity(grayimg)
    
    name = os.path.basename(image).removeprefix('out_opt_flow_registered_').removesuffix('.tif')

    if 'Atto_425' in name:       # segment Atto425 stained image
        preproc = ski.exposure.rescale_intensity(img)                                   # improves contrast

        height, width = img.shape                                                       # create contour plot
        fig, ax = plt.subplots(figsize=(width/25, height/25), dpi=25, frameon=False)    # change numbers to adjust resolution  # all numbers must be the same
        qcs = ax.contour(preproc, origin='image', levels=2)                             # qcs contains contour info
        ax.set_axis_off()
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.tight_layout(pad=0)
        fig.canvas.draw()
        plt.close(fig)

        contimg = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(height, width, 3)  # turn contour plot into image array  ## maybe try and use fig.canvas.frombuffer, as .tostring might not work for long
        contimg = ski.color.rgb2gray(contimg)
        contimg = ski.morphology.opening(contimg)

    elif 'merge' in image:        # segment merge image
        img = ski.exposure.rescale_intensity(img)

        height, width = img.shape                                                       # create contour plot
        fig, ax = plt.subplots(figsize=(width/15, height/15), dpi=15, frameon=False)    # change numbers to adjust resolution  # all numbers must be the same
        qcs = ax.contour(img, origin='image', levels=10)                                # qcs contains contour info
        ax.set_axis_off()
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.tight_layout(pad=0)
        fig.canvas.draw()
        plt.close(fig)

        contimg = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(height, width, 3)  # turn contour plot into image array  ## maybe try and use fig.canvas.frombuffer, as .tostring might not work for long
        contimg = ski.color.rgb2gray(contimg)
        contimg = ski.morphology.opening(contimg)
    
    else:                       # segment all other images
        preproc = ski.exposure.rescale_intensity(img)

        height, width = img.shape                                                       # create contour plot
        fig, ax = plt.subplots(figsize=(width/15, height/15), dpi=15, frameon=False)    # change numbers to adjust resolution  # all numbers must be the same
        qcs = ax.contour(preproc, origin='image', levels=7)                             # qcs contains contour info
        ax.set_axis_off()
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.tight_layout(pad=0)
        fig.canvas.draw()
        plt.close(fig)

        contimg = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(height, width, 3)  # turn contour plot into image array  ## maybe try and use fig.canvas.frombuffer, as .tostring might not work for long
        contimg = ski.color.rgb2gray(contimg)
        contimg = ski.morphology.opening(contimg)

    threshold = 0.9  # exact number does not matter                                     # segment contour image
    seg = ski.morphology.closing(contimg > threshold, ski.morphology.square(3))
    seg = np.invert(seg)
    seg = ski.morphology.closing(seg)

    distance = ndi.distance_transform_edt(seg)                                          # watershed
    coords = ski.feature.peak_local_max(distance, footprint=np.ones((4, 4)), labels=seg, min_distance=7)
    mask = np.zeros(img.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    wsh = ski.segmentation.watershed(-distance, markers, mask=seg)

    io.imsave(labelspath + 'labels_' + name + '.tif', wsh, check_contrast=False)


    _, num = ski.measure.label(wsh, return_num=True)                                    # label segmented regions

    # draw coordinates in plot
    plt.title(name + '   |   points of expression found: ' + str(num))

    gc.collect()
    plt.imshow(ski.color.label2rgb(wsh, grayimg, alpha=0.25, bg_label=0, bg_color=None))
    
    plt.savefig(figpath + 'Diagnostic-Figure' + name + '.png')
    
    if showplots:  # if 'showplots' is 'True'
        plt.show()
    else:
        plt.close()

    gc.collect()




def check_expression(fov, subfolder):
## checks the areas of each label in the merge image for expression in the other images
# arguments:    'fov' is the name of the current fov that should be processed
#               'subfolder' is defined in 'bashdebug.sh' and is the name of the subfolder that results are saved to
    cycles = ['c01', 'c02', 'c03', 'c04']
    nts = ['A', 'C', 'G', 'T']
    results = []        # is filled with barcode results to return
    cycleresults = []   # temp
    
    # define merged image to use as mask
    ref = io.imread('../results/Task-15/' + subfolder + '/labels/labels_merge_' + str(fov) + '.tif') 
    coords = ski.measure.regionprops_table(ref, properties=('label', 'centroid'))   # write down center coordinates of labels for returned df
    labels = np.unique(ref)                                                         # get all unique labels in 'ref'
    labels = labels[labels != 0]                                                    # exclude background

    for i in range(0, len(cycles)):                                                 # iterate over cycles
        imgs = [io.imread('../results/Task-15/' + subfolder + '/labels/labels_' + str(fov) + '_' + cycles[i] + '_Atto_425.tif'),  # define sequencing channels for this cycle
                io.imread('../results/Task-15/' + subfolder + '/labels/labels_' + str(fov) + '_' + cycles[i] + '_Alexa_488.tif'),
                io.imread('../results/Task-15/' + subfolder + '/labels/labels_' + str(fov) + '_' + cycles[i] + '_Alexa_568.tif'),
                io.imread('../results/Task-15/' + subfolder + '/labels/labels_' + str(fov) + '_' + cycles[i] + '_Alexa_647.tif')]
        check = dict(zip(nts, imgs))

        # iterate over 'labels'
        for label in labels:
            labelmask = (ref == label)                                              # create a mask for the current label

            j = 0  # flag if any nucleotides found
            k = 0  # flag if multiple nucleotides found
            for l in range(0, len(check)):                                          # iterate over images in 'check'
                if np.sum(check[nts[l]][labelmask])>0 and j == 0:                   # is there something in the mask?
                    cycleresults.append(nts[l])                                     # write corresponding NT into results
                    j = 1
                if np.sum(check[nts[l]][labelmask])>0 and j == 1 and nts[l] not in cycleresults: # if another NT was already found
                    lastentry = cycleresults[-1] + nts[l]                           # add the current NT to the already found NT
                    cycleresults[-1] = lastentry
                    k = 1
            if j == 0:                                                              # if no NTs were found
                cycleresults.append('0')                                            # 0 indicates no NT found
            if k == 1:
                lastentry = cycleresults[-1]
                cycleresults[-1] = '[' + lastentry + ']'
        
        results.append(cycleresults)
        cycleresults = []
    
    # return results for this fov
    return results, coords




#def check_expression_v2(fov, subfolder):  ## tried out, but does not perform as well as 'check_expression()'
### checks the areas of each label in the merge image for intensity peaks in the sequencing image
## arguments:    'fov' is the name of the current fov that should be processed
##               'subfolder' is defined in 'bashdebug.sh' and is the name of the subfolder that results are saved to
#    cycles = ['c01', 'c02', 'c03', 'c04']
#    nts = ['A', 'C', 'G', 'T']
#    results = []        # is filled with barcode results to return
#    cycleresults = []   # temp
#    
#    # define merged image to use as mask
#    ref = io.imread('../results/Task-15/' + subfolder + '/labels/labels_merge' + str(fov) + '.tif') 
#    coords = ski.measure.regionprops_table(ref, properties=('label', 'centroid'))  # write down center coordinates of labels for returned df
#    labels = np.unique(ref)                                                     # get all unique labels in 'ref'
#    labels = labels[labels != 0]                                                # exclude background
#
#    for i in range(0, len(cycles)):                                             # iterate over cycles
#        imgs = [io.imread('../data/selected-tiles/selected-tiles/out_opt_flow_registered' + str(fov) + cycles[i] + '_Atto_425.tif'),  # define sequencing channels for this cycle
#                io.imread('../data/selected-tiles/selected-tiles/out_opt_flow_registered' + str(fov) + cycles[i] + '_Alexa_488.tif'),
#                io.imread('../data/selected-tiles/selected-tiles/out_opt_flow_registered' + str(fov) + cycles[i] + '_Alexa_568.tif'),
#                io.imread('../data/selected-tiles/selected-tiles/out_opt_flow_registered' + str(fov) + cycles[i] + '_Alexa_647.tif')]
#        check = dict(zip(nts, imgs))
#
#        intthresh = []
#        for j in range(0, len(imgs)):                                           # set a threshold for each image based on its average intensity
#            intthresh.append(np.mean(imgs[j])*4)
#        intthresh = dict(zip(nts, intthresh))
#
#        # iterate over 'labels'
#        for label in labels:
#            labelmask = (ref == label)                                          # create a mask for the current label
#
#            j = 0  # flag if any nucleotides found
#            k = 0  # flag if multiple nucleotides found
#            for l in range(0, len(check)):                                      # iterate over images in 'check'
#                if np.mean(check[nts[l]][labelmask])>intthresh[nts[l]] and j == 0:  # is there something with higher than average intensity in the centroid?
#                    cycleresults.append(nts[l])                                 # write corresponding NT into results
#                    j = 1
#                if np.mean(check[nts[l]][labelmask])>intthresh[nts[l]] and j == 1 and nts[l] not in cycleresults: # if another NT was already found
#                    lastentry = cycleresults[-1] + nts[l]                       # add the current NT to the already found NT
#                    cycleresults[-1] = lastentry
#                    k = 1
#            if j == 0:                                                          # if no NTs were found
#                cycleresults.append('0')                                        # 0 indicates no NT found
#            if k == 1:
#                lastentry = cycleresults[-1]
#                cycleresults[-1] = '[' + lastentry + ']'
#        
#        results.append(cycleresults)
#        cycleresults = []
#    
#    # return results for this fov
#    return results, coords




def decode_barcode(fov, subfolder):
## segments all images in one fov unsing 'count_nuclei_v4()', decodes barcode using 'check_expression()'
# arguments:    'fov' is the name of the fov that should be
#               'subfolder' is the name of the subfolder that results are saved to

    # create folders
    mergepath = '../results/Task-15/' + subfolder + '/merge/'                   # merged images are saved here
    if not os.path.exists(mergepath):
            os.makedirs(mergepath)
    
    labelspath = '../results/Task-15/' + subfolder + '/labels/'                 # labels are saved here
    if not os.path.exists(labelspath):
            os.makedirs(labelspath)

    figpath = '../results/Task-15/' + subfolder + '/diagnostic-figures/'        # diagnostic figures are saved here
    if not os.path.exists(figpath):
            os.makedirs(figpath)
    
    resultspath = '../results/Task-15/' + subfolder + '/individual-results/'    # results dataframe is saved here
    if not os.path.exists(resultspath):
            os.makedirs(resultspath)
    
    # create merged image
    cycles = ['c01', 'c02', 'c03', 'c04']
    images = [] # is filled with filepaths to images specified by fovs in 'fov'
    temp = []
    imgs = []   # temp
    
    searchpattern = '../data/selected-tiles/selected-tiles/out_opt_flow_registered_' + str(fov) + '_*.tif'
    for file in glob.glob(searchpattern):
        excludepattern = ['DAPI', 'Atto_490LS']
        if excludepattern[0] not in file and excludepattern[1] not in file:
            images.append(file)

    for i in range(0, len(cycles)):
        searchpattern = str(cycles[i])
        for file in images:
            if searchpattern in file:
                imgs.append(io.imread(file, dtype=np.float64))  # read all images that should be merged
        if i == 0:  # only in cycle c01
            #merge = np.sum(imgs, axis=0)
            merge = np.max(imgs, axis=0)
    
    # merge images
    #merge = np.sum(imgs, axis=0)                                
    temp = mergepath + 'merge_' + str(fov) + '.tif'
    io.imsave(temp, merge, check_contrast=False)
    images.append(temp)

    gc.collect()  # delete unused variables

    # segment nuclei and decode barcode using decode_barcode()
    results = []    # is filled with results returned from 'check_expression()'

    # segment all images in this fov
    for i in range(0, len(images)):
        count_nuclei_v4(images[i], labelspath=labelspath, figpath=figpath, showplots=False)
        gc.collect()
    
    # measure intensities, decode barcodes and create results dataframe
    results, coords = check_expression(fov, subfolder=subfolder)  # get NTs

    barcodes = [''.join(barcode) for barcode in zip(*results)]  # reformat barcodes for df

    tempdf = pd.DataFrame(columns=['fov'])

    taglist = pd.read_csv('../data/taglist.csv')                # get taglist that contains genes and barcodes
    genes = []
    for barcode in barcodes:
        if '0' in barcode or '[' in barcode:                    # exclude invalid barcodes
            genes.append('invalid')
        else:                                                   # check all valid barcodes
            j = 0 # flag if gene found
            for i in range(0,len(taglist.index)):               # iterate over rows in taglist
                if taglist.at[i, 'Code'] == barcode:
                    genes.append(taglist.at[i, 'Name'])         # append corresponding gene to 'gene'
                    j = 1 
            if j == 0:                                          # if no genes found
                genes.append('no match')

    fovname = fov.strip('_')
    totemp = dict(zip(['x', 'y', 'barcode', 'gene'], [coords['centroid-0'], coords['centroid-1'], barcodes, genes]))
    tempdf2 = pd.DataFrame(data=totemp)
    tempdf = pd.concat([tempdf, tempdf2])
    tempdf['fov'] = tempdf['fov'].fillna(str(fovname))
    tempdf.to_csv(resultspath + 'results_' + str(fovname) + '.csv')

    gc.collect()  # delete unused variables

    return tempdf




def count_genes(subfolder, omit=False, saveconcat=True):
## reads the .csv file resulting from 'decode_barcodes()' and counts identified genes
# arguments:    'subfolder' is the name of the subfolder that contains the results from 'decode_barcodes()'
#               'omit=True' omits all 'invalid' and 'no match' entries from the returned dataframe
#               'saveconcat=True' saves the concatenated dataframe that contains the results from 'decode_barcodes()'
#                   if this is 'False', only the gene count dataframe is saved

    # concatenate all results files from 'decode_barcodes()'
    tempresults = []                                        # get all results files from 'decode_barcodes()'
    searchpattern = '../results/Task-15/' + subfolder + '/individual-results/results_*.csv'
    for file in glob.glob(searchpattern):
        tempresults.append(pd.read_csv(file))
    concatdf = pd.concat(tempresults, ignore_index=True)    # concatenated df
    del concatdf['Unnamed: 0']

    if saveconcat:  # if 'saveconcat' is 'True
        concatdf.to_csv('../results/Task-15/' + subfolder + '/concatenated-results.csv')

    columnnames = ['fov']  # other columns are added later
    gcdf = pd.DataFrame(columns=columnnames)                # gene count df to be returned
    decode_byfov = concatdf.groupby(concatdf['fov'])

    for name in concatdf['fov'].unique():                   # fill 'df' with data
        temp = decode_byfov.get_group(name)
        counts = temp['gene'].value_counts()
        gcdf = pd.concat([gcdf, counts])
        gcdf = gcdf.fillna(value=name)
    
    if omit:                                                # remove 'invalid' and 'no match' from df
        gcdf = gcdf[gcdf.index != 'invalid']
        gcdf = gcdf[gcdf.index != 'no match']

    gcdf = gcdf.rename(columns={0: 'count'})                # make df look nice
    gcdf = gcdf.reset_index(names='gene')
    gcdf = gcdf[['fov', 'gene', 'count']]

    gcdf.to_csv('../results/Task-15/' + subfolder + '/gene-counts.csv')

    return gcdf




def visualize(fov, subfolder, markers='shapes', labelonlygenes=True):
## launches Napari viewer, adds all channels and all cycles of one fov and marks points with identified barcode
# arguments:    'fov' is the name of the fov that should be displayed
#               'subfolder' is the name of the subfolder where the results tables are
#               'markers' decides how the identified points should be marked
#                   'markers='shapes'' marks the points with polygon outlines
#                   'markers='points'' marks the points with points at their centroid position
#               'labelonlygenes=True' only displays text labels for identified genes and omits labels for 'invalid' and 'no match' results
    cycles      = ['c01', 'c02', 'c03', 'c04']
    channels    = ['C', 'G', 'T', 'A', 'DAPI']
    colors      = ['red', 'green', 'yellow', 'magenta', 'blue']
    opacities   = [1, 1, 1, 1, 0.5]
    opacities   = opacities*len(cycles)
    imgs = []

    # get all images
    searchpattern = '../data/selected-tiles/selected-tiles/out_opt_flow_registered_' + fov + '_*.tif'
    excludepattern = 'Atto_490LS'
    for file in glob.glob(searchpattern):
        if excludepattern not in file:
            #imgs.append(file)
            imgs.append(io.imread(file))

    # sort images into dict and create a multidimensional image
    height, width = imgs[0].shape[:2]                   # get shape of the images
    displaylayers = np.zeros((4, 5, height, width), dtype=imgs[0].dtype)
    i = 0
    for j in range(0, len(cycles)):
        for k in range(0, len(channels)):
            displaylayers[j, k] = imgs[i]
            i += 1

    barcodes = pd.read_csv('../results/Task-15/' + subfolder + '/individual-results/results_' + fov + '.csv', usecols=['barcode', 'gene'])

    if labelonlygenes:
        textlabels = []
        for i in range(0, barcodes['barcode'].count()): # display text only for identified genes
            if barcodes['gene'][i] != 'invalid' and barcodes['gene'][i] != 'no match':
                textlabels.append(barcodes['barcode'][i] + ':\n' + barcodes['gene'][i])
            else: 
                textlabels.append('')
        features = {'textlabels': textlabels}
        text = {
            'string': '{textlabels}',
            'size': 5,
            'color': 'white',
            'translation': np.array([-2, 0]),
            'anchor': 'upper_left'
        }
    if not labelonlygenes:
        features = {'barcode': barcodes['barcode'], 'gene': barcodes['gene']}
        text = {
            'string': '{barcode}:\n{gene}',
            'size': 5,
            'color': 'white',
            'translation': np.array([-2, 0]),
            'anchor': 'upper_left'
        }

    # launch and fill napari viewer
    viewer = napari.view_image(                         # launch and add multidimensional image
        data=displaylayers, 
        name=channels, 
        blending='additive', 
        channel_axis=1, 
        colormap=colors,
        opacity=opacities)

    if markers=='shapes':
        # get contours to annotate each point of expression
        labeledmask = io.imread('../results/Task-15/' + subfolder + '/labels/labels_merge_' + fov + '.tif')
        poes = np.unique(labeledmask)  # get unique labels as poes
        poes = poes[poes != 0]
        contours = []
        for poe in poes:
            mask = labeledmask == poe
            mask = ski.morphology.remove_small_holes(mask, area_threshold=500)  # fill holes to ensure only one contour per poe, otherwise text labeling doesn't work
            contours.extend(ski.measure.find_contours(mask))
        
        # add to napari
        viewer.add_shapes(                              # add outline shapes
            contours, 
            shape_type='polygon', 
            edge_color='white', 
            face_color='transparent', 
            name='Outlines',
            opacity=0.5)
        viewer.add_shapes(                              # add outline text
            contours, 
            shape_type='polygon', 
            edge_color='transparent', 
            face_color='transparent', 
            name='Outlines-Text', 
            features=features, 
            text=text, 
            opacity=1)
    
    if markers == 'points':
        # get coorinates for points to annotate each point of expression
        points = pd.read_csv('../results/Task-15/' + subfolder + '/individual-results/results_' + fov + '.csv', usecols=['x', 'y'])

        # add to napari
        viewer.add_points(points, size=2, border_width=0, features=features, text=text)

    return viewer