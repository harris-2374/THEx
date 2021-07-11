
def remove_heterotachy_info(l):
    if "[" not in l:
        return l
    else:
        open_brackets = [i for i, x in enumerate(l) if x == "["]
        close_brackets = [i for i, x in enumerate(l) if x == "]"]
        final_string = f'{l[:open_brackets[0]]}'
        for ob, cb in zip(open_brackets[1:], close_brackets[:-1]):
            final_string += l[cb+1:ob]
        final_string += l[close_brackets[-1]+1:]
        return final_string

def get_first_tree(f):
    with open(f) as fh:
        for l in fh:
            return l.strip()

def clean_iqtree_output():
    treefile_dir = path("")
    chromdirs = [d for d in treefile_dir.iterdir() if d.is_dir()]
    
    # Get first line (tree) from file
    for d in chromdirs:
        files = [f for f in d.iterdir() if f.is_file()]
        for f in files:
            newickTree = get_first_tree(f)
            # Remove heterotacky info
            cleanTree = remove_heterotachy_info
            print(cleanTree)
            # Overwrite original file
            # with open(tf, 'w') as fh:
            #     fh.write(final_string)
            #     continue
    return

if __name__ == '__name__':
    clean_iqtree_output()