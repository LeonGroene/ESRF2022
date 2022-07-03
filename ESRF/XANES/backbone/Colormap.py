"""

:platform: Unix, Windows 
:last changed: 2022-05-13

.. moduleauthor:: Tjark Leon Raphael Gröne <tgroene@physnet.uni-hamburg.de>


"""
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib import cm as cmm

#################################### Color Map #######################################

def cm(N):
    # light Bue --> White
    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(62/256, 1, N)
    vals2[:, 1] = np.linspace(204/256, 1, N)
    vals2[:, 2] = np.linspace(0, 0, N)
    cmap1 = ListedColormap(vals2)

    # light Yellow --> Yellow
    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(0, 62/256, N)
    vals2[:, 1] = np.linspace(1, 204/256, N)
    vals2[:, 2] = np.linspace(238/256, 0, N)
    cmap2 = ListedColormap(vals2)

    # Yellow -->  Orange
    vals3 = np.ones((N, 4))
    vals3[:, 0] = np.linspace(0, 0, N)
    vals3[:, 1] = np.linspace(164/256, 1, N)
    vals3[:, 2] = np.linspace(1, 238/256, N)
    cmap3 = ListedColormap(vals3)

    # Orange --> Red
    vals4 = np.ones((N, 4))
    vals4[:, 0] = np.linspace(58/256, 0, N)
    vals4[:, 1] = np.linspace(49/256, 164/256, N)
    vals4[:, 2] = np.linspace(111/256, 1, N)
    cmap4 = ListedColormap(vals4)

    # Blue --> light Blue 
    vals = np.ones((N, 4))
    vals[:, 0] = np.linspace(0, 58/256, N)
    vals[:, 1] = np.linspace(0, 49/256, N)
    vals[:, 2] = np.linspace(0, 111/256, N)
    cmap5 = ListedColormap(vals)



    # Stack Blue --> Red
    newcolors = np.vstack((
                            # cmap5(np.linspace(0, 1, 128)),
                            cmap4(np.linspace(0, 1, int(N/4))), 
                            cmap3(np.linspace(0, 1, int(N/4))),
                            cmap2(np.linspace(0, 1, int(N/4))),
                            cmap1(np.linspace(0, 1, int(N/4)))
                            ))
    

    colormap = []
    for i in range(0,len(newcolors)):
        colormap.append([newcolors[i][0],newcolors[i][1],newcolors[i][2]])
    
    return colormap


def cmRB_Time(N):
    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(0, 0.55, N)
    vals2[:, 1] = np.linspace(0, 1, N)
    vals2[:, 2] = np.linspace(1, 1, N)
    cmap1 = ListedColormap(vals2)

    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(0.55, 1, N)
    vals2[:, 1] = np.linspace(1, 0, N)
    vals2[:, 2] = np.linspace(1, 0, N)
    cmap2 = ListedColormap(vals2)

    newcolors = np.vstack((
                        cmap1(np.linspace(0, 1, int((3*N)/8))),
                        cmap2(np.linspace(0, 1, int((5*N)/8))),
                        ))
    colormap = []
    for i in range(0,len(newcolors)):
        colormap.append([newcolors[i][0],newcolors[i][1],newcolors[i][2]])
    
    return ListedColormap(colormap, name='RedBlue')

def cmRB(N):
    return cmm.get_cmap('bwr', N)

linestyle_tuple = [
     (0, (1, 10)),
     (0, (1, 1)),
     (0, (1, 1)),
     (0, (5, 10)),
     (0, (5, 5)),
     (0, (5, 1)),
     (0, (3, 10, 1, 10)),
     (0, (3, 5, 1, 5)),
     (0, (3, 1, 1, 1)),
     (0, (3, 5, 1, 5, 1, 5)),
     (0, (3, 10, 1, 10, 1, 10)),
     (0, (3, 1, 1, 1, 1, 1)),
     (0, (1, 10)),
     (0, (1, 1)),
     (0, (1, 1)),
     (0, (5, 10)),
     (0, (5, 5)),
     (0, (5, 1)),
     (0, (3, 10, 1, 10)),
     (0, (3, 5, 1, 5)),
     (0, (3, 1, 1, 1)),
     (0, (3, 5, 1, 5, 1, 5)),
     (0, (3, 10, 1, 10, 1, 10)),
     (0, (3, 1, 1, 1, 1, 1)),]

    
color_order = [ "black", "red", "blue", "lawngreen", 
                "fuchsia", "grey", "yellow","brown",
                "cyan", "orchid", "olive","deeppink",
                "mediumspringgreen", "goldenrod","steelblue","firebrick",
                "forestgreen","indigo", "peru", "navy",
                "dimgray", "darkorange", "green", "orangered",
                "slategray", "orchid", "tomato", "aquamarine",
                "tan", "maroon", "plum", "greenyellow"]		

def cm2d():
    # RGB Colorrange
    N = 256

    
    # light Bue --> White
    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(62/256, 1, N)
    vals2[:, 1] = np.linspace(204/256, 1, N)
    vals2[:, 2] = np.linspace(0, 0, N)
    cmap1 = ListedColormap(vals2)

    # light Yellow --> Yellow
    vals2 = np.ones((N, 4))
    vals2[:, 0] = np.linspace(0, 62/256, N)
    vals2[:, 1] = np.linspace(1, 204/256, N)
    vals2[:, 2] = np.linspace(238/256, 0, N)
    cmap2 = ListedColormap(vals2)

    # Yellow -->  Orange
    vals3 = np.ones((N, 4))
    vals3[:, 0] = np.linspace(0, 0, N)
    vals3[:, 1] = np.linspace(164/256, 1, N)
    vals3[:, 2] = np.linspace(1, 238/256, N)
    cmap3 = ListedColormap(vals3)

    # Orange --> Red
    vals4 = np.ones((N, 4))
    vals4[:, 0] = np.linspace(58/256, 0, N)
    vals4[:, 1] = np.linspace(49/256, 164/256, N)
    vals4[:, 2] = np.linspace(111/256, 1, N)
    cmap4 = ListedColormap(vals4)

    # Blue --> light Blue 
    vals = np.ones((N, 4))
    vals[:, 0] = np.linspace(0, 58/256, N)
    vals[:, 1] = np.linspace(0, 49/256, N)
    vals[:, 2] = np.linspace(0, 111/256, N)
    cmap5 = ListedColormap(vals)



    # Stack Blue --> Red
    newcolors = np.vstack((
                            cmap5(np.linspace(0, 1, 128)),
                            cmap4(np.linspace(0, 1, 128)), 
                            cmap3(np.linspace(0, 1, 128)),
                            cmap2(np.linspace(0, 1, 128)),
                            cmap1(np.linspace(0, 1, 512))
                            ))
    newcmp = ListedColormap(newcolors, name='RedBlue')
    return newcmp