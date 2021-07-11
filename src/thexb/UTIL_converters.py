def convert_window_size_to_int(ws):
    """This function converts the shorthand input window size
    and returns an integer of the same value (i.e. "100kb" == int(100000))"""
    ws = ws.lower()
    window_units = ws[-2:]
    try:
        window_int = int(ws[:-2])
    except ValueError:
        raise ValueError("Please give window size as integer (i.e., 10kb)")
    
    converters = {
        'bp': int(ws[:-2]),
        'kb': int(ws[:-2])*(10**3),
        'mb': int(ws[:-2])*(10**6),
        'gb': int(ws[:-2])*(10**9)   
    }

    try:
        return converters[window_units]
    except KeyError:
        print(window_units)
        raise KeyError("Invalid unit type provided. Options=[bp, kb, mb, gb]")
        # raise TypeError("Invalid window size value, must be a integer followed by bp, kb, mb, gb (i.e., 10bp, 10kb, 10mb, 10gb)")
    
