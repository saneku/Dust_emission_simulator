from matplotlib.colors import ListedColormap

def make_colormap(cmap_arr, name):
    """
    Make a matplotlib colormap based on a ncview-type definition of a colormap

    source: https://stackoverflow.com/a/19561196
    """
    ctmp = []
    i = 0
    while i < len(cmap_arr):
        ctmp.append([cmap_arr[i]/255., cmap_arr[i+1]/255., cmap_arr[i+2]/255.])
        i += 3

    return ListedColormap(ctmp, name=name, N=None)
