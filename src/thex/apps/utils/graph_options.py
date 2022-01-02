def color_swatches():
    color_swatches ={
        'Plotly': [
            '#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52'
            ]
        , 
        'D3': [
            '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF'
            ]
        , 
        'G10': [
            '#3366CC', '#DC3912', '#FF9900', '#109618', '#990099', '#0099C6', '#DD4477', '#66AA00', '#B82E2E', '#316395'
            ]
        , 
        'T10': [
            '#4C78A8', '#F58518', '#E45756', '#72B7B2', '#54A24B', '#EECA3B', '#B279A2', '#FF9DA6', '#9D755D', '#BAB0AC'
            ]
        , 
        'Alphabet': [
            '#AA0DFE', '#3283FE', '#85660D', '#782AB6', '#565656', '#1C8356', '#16FF32', '#F7E1A0', '#E2E2E2', '#1CBE4F', 
            '#C4451C', '#DEA0FD', '#FE00FA', '#325A9B', '#FEAF16', '#F8A19F', '#90AD1C', '#F6222E', '#1CFFCE', '#2ED9FF', 
            '#B10DA1', '#C075A6', '#FC1CBF', '#B00068', '#FBE426', '#FA0087'
            ]
        , 
        'Dark24': [
            '#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', 
            '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', 
            '#DA60CA', '#6C4516', '#0D2A63', '#AF0038'
            ]
        , 
        'Light24': [
            '#FD3216', '#00FE35', '#6A76FC', '#FED4C4', '#FE00CE', '#0DF9FF', '#F6F926', '#FF9616', '#479B55', '#EEA6FB', 
            '#DC587D', '#D626FF', '#6E899C', '#00B5F7', '#B68E00', '#C9FBE5', '#FF0092', '#22FFA7', '#E3EE9E', '#86CE00', 
            '#BC7196', '#7E7DCD', '#FC6955', '#E48F72'
            ]
        , 
        'Set1': [
            'rgb(228,26,28)', 'rgb(55,126,184)', 'rgb(77,175,74)', 'rgb(152,78,163)', 'rgb(255,127,0)', 'rgb(255,255,51)', 
            'rgb(166,86,40)', 'rgb(247,129,191)', 'rgb(153,153,153)'
            ]
        , 
        'Pastel1': [
            'rgb(251,180,174)', 'rgb(179,205,227)', 'rgb(204,235,197)', 'rgb(222,203,228)', 'rgb(254,217,166)', 'rgb(255,255,204)', 
            'rgb(229,216,189)', 'rgb(253,218,236)', 'rgb(242,242,242)'
            ]
        , 
        'Dark2': [
            'rgb(27,158,119)', 'rgb(217,95,2)', 'rgb(117,112,179)', 'rgb(231,41,138)', 'rgb(102,166,30)', 'rgb(230,171,2)', 
            'rgb(166,118,29)', 'rgb(102,102,102)'
            ]
        , 
        'Set2': [
            'rgb(102,194,165)', 'rgb(252,141,98)', 'rgb(141,160,203)', 'rgb(231,138,195)', 'rgb(166,216,84)', 'rgb(255,217,47)', 
            'rgb(229,196,148)', 'rgb(179,179,179)'
        ], 
        'Pastel2': [
            'rgb(179,226,205)', 'rgb(253,205,172)', 'rgb(203,213,232)', 'rgb(244,202,228)', 'rgb(230,245,201)', 'rgb(255,242,174)', 
            'rgb(241,226,204)', 'rgb(204,204,204)'
        ], 
        'Set3': [
            'rgb(141,211,199)', 'rgb(255,255,179)', 'rgb(190,186,218)', 'rgb(251,128,114)', 'rgb(128,177,211)', 'rgb(253,180,98)', 
            'rgb(179,222,105)', 'rgb(252,205,229)', 'rgb(217,217,217)', 'rgb(188,128,189)', 'rgb(204,235,197)', 'rgb(255,237,111)'
        ], 
        'Antique': [
            'rgb(133, 92, 117)', 'rgb(217, 175, 107)', 'rgb(175, 100, 88)', 'rgb(115, 111, 76)', 'rgb(82, 106, 131)', 'rgb(98, 83, 119)', 
            'rgb(104, 133, 92)', 'rgb(156, 156, 94)', 'rgb(160, 97, 119)', 'rgb(140, 120, 93)', 'rgb(124, 124, 124)'
        ], 
        'Bold': [
            'rgb(127, 60, 141)', 'rgb(17, 165, 121)', 'rgb(57, 105, 172)', 'rgb(242, 183, 1)', 'rgb(231, 63, 116)', 'rgb(128, 186, 90)', 
            'rgb(230, 131, 16)', 'rgb(0, 134, 149)', 'rgb(207, 28, 144)', 'rgb(249, 123, 114)', 'rgb(165, 170, 153)'
        ], 
        'Pastel': [
            'rgb(102, 197, 204)', 'rgb(246, 207, 113)', 'rgb(248, 156, 116)', 'rgb(220, 176, 242)', 'rgb(135, 197, 95)', 'rgb(158, 185, 243)', 
            'rgb(254, 136, 177)', 'rgb(201, 219, 116)', 'rgb(139, 224, 164)', 'rgb(180, 151, 231)', 'rgb(179, 179, 179)'
        ], 
        'Prism': [
            'rgb(95, 70, 144)', 'rgb(29, 105, 150)', 'rgb(56, 166, 165)', 'rgb(15, 133, 84)', 'rgb(115, 175, 72)', 'rgb(237, 173, 8)', 
            'rgb(225, 124, 5)', 'rgb(204, 80, 62)', 'rgb(148, 52, 110)', 'rgb(111, 64, 112)', 'rgb(102, 102, 102)'
        ], 
        'Safe': [
            'rgb(136, 204, 238)', 'rgb(204, 102, 119)', 'rgb(221, 204, 119)', 'rgb(17, 119, 51)', 'rgb(51, 34, 136)', 'rgb(170, 68, 153)', 
            'rgb(68, 170, 153)', 'rgb(153, 153, 51)', 'rgb(136, 34, 85)', 'rgb(102, 17, 0)', 'rgb(136, 136, 136)'
        ], 
        'Vivid': [
            'rgb(229, 134, 6)', 'rgb(93, 105, 177)', 'rgb(82, 188, 163)', 'rgb(153, 201, 69)', 'rgb(204, 97, 176)', 'rgb(36, 121, 108)', 
            'rgb(218, 165, 27)', 'rgb(47, 138, 196)', 'rgb(118, 78, 159)', 'rgb(237, 100, 90)', 'rgb(165, 170, 153)'
        ]}
    return color_swatches


def figure_output_scales():
    return [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]


def graph_templates():
    return ['plotly', 'plotly_white', 'plotly_dark', 'ggplot2', 'presentation', 'seaborn', 'simple_white']


def snapshot_file_type():
    return ['svg', 'jpeg', 'png', 'webp']


def marker_width_options():
    options = ['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.0', '2.25', '2.5']
    return [{'label': s, 'value': s} for s in options]


def font_size_options():
    return [{'label': s, 'value': s} for s in range(5, 50, 1)]

def font_families():
    return ["Arial", "Balto", "Courier New", 
            "Droid Sans", "Droid Serif", "Droid Sans Mono", 
            "Gravitas One", "Old Standard TT", "Open Sans", 
            "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]