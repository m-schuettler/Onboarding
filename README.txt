README 
2024-09-02 - Magdalena Schüttler

Onboarding Base Module

DATA

The dataset consists of 171 selected tiles from the right hemisphere of a mouse brain. Each tile has the dimensions 1000x1000 pixels and contains 6 imaging channels (nuclei channel, anchor channel and 4 coding channels) in which different markers were used. For each tile 4 sequencing rounds were performed.

The dataset consists of the following files:
'channel_infos.csv'         - Relates each image channel to its corresponding code
'taglist.csv'               - Codebook of barcodes used in the experiment
'tile_names.csv'            - Names (coordinates) of selected tiles
'selected_tiles_map.png'    - Map of selected tiles
'selected-tiles.zip'        - Registered tif images of selected tiles
    - X.._Y..                   - Field of View: position of the image as part of the entire brain hemisphere image
    - c..                       - cycle
    - Alexa_../Atto_../DAPI     - dye
        -- DAPI                     - nucleus
        -- Atto_490LS               - anchor
        -- Atto_425                 - A (Adenosine)
        -- Alexa_488                - C (Cytosine)
        -- Alexa_568                - G (Guanosine)
        -- Alexa_647                - T (Thymidine)
'decoding.zip'              - Decoding of the selected ISS tiles via different methods