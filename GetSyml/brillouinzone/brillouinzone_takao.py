from __future__ import print_function
from builtins import range
import numpy as np
from collections import defaultdict


#http://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()
    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)
    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])
    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])


def get_missing_point(tr, p1, p2):
    """
    tr is a list of 3 indices (the indices of the vertices of a given triangle).
    p1 and p2 must be in the tr list, the third index is returned.

    ValueError is raised if more than one entry is not among p1 and p2,
    or if all points are among p1 and p2 (e.g. a point repeated twice)
    """
    missing = None
    for idx, vertex in enumerate(tr):
        if vertex == p1 or vertex == p2:
            pass
        else:
            if missing is not None:
                raise ValueError("Two missing points found...")
            missing = vertex

    if missing is None:
        raise ValueError("No missing points!")

    return missing

def are_coplanar(v1,v2,v3):
    """
    v1, v2, v3: 3D vectors. 
    Return True if they are coplanar, False otherwise
    """
    return float(abs(np.dot(np.cross(v1,v2),v3))) < 1.e-6


def get_BZ(b1,b2,b3):
    """
    Get the faces of the BZ given the three vectors b1,b2,b3.
    Return both triangular faces (oriented) and flat faces (if your
        plotting library prefers these - these are not oriented for the
        time being)
    """
    from scipy.spatial import Voronoi, ConvexHull, Delaunay
    ret_data = {}

    supercell_size = 3 # Is this enough?

    # G-vectors
    points3d = []
    central_idx = None
    for i in range(-supercell_size, supercell_size+1):
        for j in range(-supercell_size, supercell_size+1):
            for k in range(-supercell_size, supercell_size+1):
                if i==0 and j==0 and k==0:
                    central_idx = len(points3d)
                points3d.append(i*np.array(b1) + j*np.array(b2) + k*np.array(b3))

    # Get Voronois
    vor3d = Voronoi(np.array(points3d))
    # Get the vertices of the central" Voronoi( around the origin G=0)
    central_voronoi_3d = np.array([vor3d.vertices[idx] for idx in vor3d.regions[vor3d.point_region[central_idx]]])

    # Get the convex hull of these points (all triangular faces)
    hull = ConvexHull(central_voronoi_3d)

    #print('hull=',hull)
    ## REORIENT TRIANGLES
    ## NOTE! TODO: I should do the same for faces
    ret_data['triangles_vertices'] = hull.points.tolist()


# ################################################
#     print('lll=', len(hull.points))
#     mmm=np.matrix([[0e0]*len(hull.points)]*len(hull.points))
#     print(mmm)
#     ii=0
#     for i in hull.points:
#         ii=ii+1
#         jj=0
#         for j in hull.points:
#             jj=jj+1
#             mmm[ii-1,jj-1]=np.dot(i,j)
#     print(mmm)
#     sys.exit()


    ## Naive one
    # ret_data['triangles'] = hull.simplices.tolist()
    ## Instead, I orient them all
    ret_data['triangles'] = []
    for simplex in hull.simplices:
        points = np.array([hull.points[i] for i in simplex])
        center = points.sum(axis=0) / float(len(points))

        # Get normal vector (with direction!)
        normal = np.cross(center-points[1], center-points[0])
        # Normalize, then rescale to a small value
        normal = normal / np.linalg.norm(normal)  
        max_length = np.sqrt((points**2).sum(axis=1).max())
        normal /= max_length 
        normal *= 1.e-4
        point_up = center + normal
        point_down = center - normal
        delaunay = Delaunay(hull.points)
        is_up_inside = delaunay.find_simplex(point_up)>=0
        is_down_inside = delaunay.find_simplex(point_down)>=0

        if is_up_inside and not is_down_inside:
            correct_orientation = True
        elif not is_up_inside and is_down_inside:
            correct_orientation = False
        else:
            inside_outside_string = "inside" if is_up_inside else "outside"
            print ("WARNING! Both vectors are {}..."\
                " not changing orientation".format(inside_outside_string))
            correct_orientation = True

        if correct_orientation:
            ret_data['triangles'].append(simplex.tolist())
        else:
            ret_data['triangles'].append(simplex[::-1].tolist())


    #print hull.area, hull.volume

    # Get edge-sharing faces
    # edges has as key a tuple (sorted) with the indices of the two vertices of
    # the shared edge; the value are the indices of the triangles    
    edges = defaultdict(list)
    for simplex_idx, simplex in enumerate(hull.simplices):
        edges[tuple(sorted([simplex[0],simplex[1]]))].append(simplex_idx)
        edges[tuple(sorted([simplex[1],simplex[2]]))].append(simplex_idx)
        edges[tuple(sorted([simplex[2],simplex[0]]))].append(simplex_idx)
    # convert to dictionary of lists (from defaultdict of sets)
    edges = dict(edges) #{k: v for k, v in edges.iteritems()}

    ### Create now the list of faces, merging the triangles that share an
    ### edge and are coplanar. Note: this works only if up to two triangles
    ### must be merged; if three or more, this will not produce the expected
    ### result
    #print edges

   
    # Store merge operations to perform
    merge_with = defaultdict(set)

    # List of found simplices that have been merged; will be used at the 
    # end to add the triangles that have not been merged, if any
    merged_simplices = []
    for (p1, p2), triangles in edges.items():
        # I do it many times, but it's the easiest way (and anyway it's
        # a set, so it should be fast: I add a point to be merged with 
        # itself
        merge_with[triangles[0]].add(triangles[0])
        merge_with[triangles[1]].add(triangles[1])

        if len(triangles) != 2:
            # An edge shared by less (or more) than 2 triangles?
            print("Warning!", p1, p2, triangles)
            continue
        else:
            # Check if two triangles are coplanar: get the other two
            # vertices that are not on the shared edge
            otherpoint0 = get_missing_point(hull.simplices[triangles[0]], p1, p2)
            otherpoint1 = get_missing_point(hull.simplices[triangles[1]], p1, p2)
            # The actual vector coordinates
            otherpoint0_p = hull.points[otherpoint0]
            otherpoint1_p = hull.points[otherpoint1]
            p1_p = hull.points[p1]
            p2_p = hull.points[p2]

            # Check if they are coplanar
            if are_coplanar(p2_p - p1_p, otherpoint0_p - p1_p, otherpoint1_p - p1_p):
                merge_with[triangles[0]].add(triangles[1])
                merge_with[triangles[1]].add(triangles[0])


    # PROBLEM TO SOLVE: we have to put together all groups
    ## E.g. we have now:
    #0: [0, 3]
    #1: [1, 2, 3]
    #2: [1, 2]
    #3: [0, 1, 3]
    ## We should get instead for all [0,1,2,3]
    # So we have to do a pass to merge them all
    # The algorithm below is probably wrong (actually, it is only 
    # probably inefficient)

    #for k, v in merge_with.iteritems():
    #    print "{}: {}".format(k, list(v))

    # Iterate untile convergence - not sure this is correct
    has_changed = True
    while has_changed:
        has_changed = False
        # Add the missing ones
        for tr in range(len(hull.simplices)):
            for other1 in merge_with[tr]:
                for other2 in merge_with[tr]:
                    if other1 not in merge_with[other2]:
                        has_changed = True
                        merge_with[other2].add(other1)
                    if other2 not in merge_with[other1]:
                        has_changed = True
                        merge_with[other1].add(other2)

    # convert to dict, and most importantly convert to list and sort
    merge_with = {k: sorted(v) for k, v in merge_with.items()}
        
    # Assign to the smallest integer idx
    merge_group = {k: v[0] for k, v in merge_with.items()}

    groups = defaultdict(list)
    # I create a reverse index
    for k, v in merge_group.items():
        groups[v].append(k)

    # List of faces (elements are lists of vertex ids)
    faces = [] 

    for group in groups.values():
        if len(group) == 1:
            faces.append([hull.points[point_idx] 
                          for point_idx in hull.simplices[group[0]]])
        else:
            # Get all points
            all_points_idx = sorted(set(
                np.concatenate([hull.simplices[g] for g in group])))

            all_points_coords = [hull.points[point_idx]
                                 for point_idx in all_points_idx]
            # Find projection in 2D: I first choose a first vector (between
            # two points v1; I find the orthogonal vector to the plane b; 
            # I find a second vector v2 orthogonal to v1 and on the plane 
            # (<=> orthogonal to v1 and b); I normalize v1 and v2; 
            # I get the components of the vectors w.r.t. v1 and v2
            # NOTE: there is at least 1 triangle => at least 3 points
            v1 = all_points_coords[1] - all_points_coords[0]
            temp_v2 = all_points_coords[2] - all_points_coords[0]
            b = np.cross(v1, temp_v2)
            v2 = np.cross(v1, b)
            # Normalize
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            # get components
            x = [np.dot(point, v1) for point in all_points_coords]
            y = [np.dot(point, v2) for point in all_points_coords]
            # 2. do convexhull in 2D
            hull_face2d = ConvexHull(np.array([x,y]).T)
            # 3. get point indices of convex hull, convert back to original 
            # 3D indices
            # 2dhull.vertices contains the segments, but with ids in the
            # smaller subset. We want the indices in the initial set:
            actual_points_idx = [all_points_idx[subset_idx] for 
                                 subset_idx in hull_face2d.vertices]
            # 4. add to faces list
            faces.append([hull.points[point_idx].tolist()
                          for point_idx in actual_points_idx])
        
    ret_data['faces'] = faces
    return ret_data,points3d




def plotws(b1,b2,b3):
    from pylab import figure, show
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    #draw a vector
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
            FancyArrowPatch.draw(self, renderer)

    faces_data,points3d = get_BZ(b1,b2,b3)
    import json
    #print(json.dumps(faces_data))

    faces_coords = faces_data['faces']

    faces_count = defaultdict(int)
    for face in faces_coords:
        faces_count[len(face)] += 1

    for num_sides in sorted(faces_count.keys()):
        print("{} faces: {}".format(num_sides, faces_count[num_sides]))

    fig = figure() #figsize=figaspect(1))
    ax = fig.add_subplot(111, projection='3d')
    ax.add_collection3d(Poly3DCollection(faces_coords,linewidth=1, alpha=0.3, edgecolor="k", facecolor="#ccccff"))

    #draw origin
    ax.scatter([0],[0],[0],color="g",s=100)

    axes_length = 1 

# Add axes
    ax.add_artist(Arrow3D((0,axes_length),(0,0),(0,0), 
          mutation_scale=20, lw=2, arrowstyle="-|>", color="k"))
    ax.add_artist(Arrow3D((0,0),(0,axes_length),(0,0),
          mutation_scale=20, lw=2, arrowstyle="-|>", color="k"))
    ax.add_artist(Arrow3D((0,0),(0,0),(0,axes_length),
         mutation_scale=20, lw=2, arrowstyle="-|>", color="k"))
    ax.text(axes_length,0,0, 'ex', (1,0,0))
    ax.text(0,axes_length,0, 'ey', (1,0,0))
    ax.text(0,0,axes_length, 'ez', (1,0,0))

# three primitive vectors    
    ax.add_artist(Arrow3D((0,b1[0]),(0,b1[1]),(0,b1[2]), 
          mutation_scale=21, lw=1, arrowstyle="-|>", color="r"))
    ax.add_artist(Arrow3D((0,b2[0]),(0,b2[1]),(0,b2[2]),
          mutation_scale=21, lw=1, arrowstyle="-|>", color="r"))
    ax.add_artist(Arrow3D((0,b3[0]),(0,b3[1]),(0,b3[2]),
          mutation_scale=21, lw=1, arrowstyle="-|>", color="r"))
    ax.text(b1[0], b1[1], b1[2], 'qlat1', (1,0,0))
    ax.text(b2[0], b2[1], b2[2], 'qlat2', (1,0,0))
    ax.text(b3[0], b3[1], b3[2], 'qlat3', (1,0,0))
##

    import sys,re
    iii=''
    eee=''
    x=['']*3
    y=['']*3
    n=''

    symlfile='syml.'+sys.argv[1]
    sfile = open(symlfile,'r').read().split('\n')

# look for nline
#    print ('sfile=',sfile)
    i=0
    for iline in sfile:
        i=i+1
        ilr=re.split('\s+',iline)
        n=int(ilr[0])
        if(n==0): break
    nline=i-1
    #print( 'nline=',nline)
    for iline in sfile:
        i=i+1
        #print(i,iline)
        ilr=re.split('\s+',iline)
        n=int(ilr[0])
        if(n==0): break
        x= [float(ilr[i]) for i in range(1,4)]
        y= [float(ilr[i]) for i in range(4,7)]
        iii = ilr[7]
        eee = ilr[8]
        #print( n,x[0:3],y[0:3],iii,eee)
        ax.add_artist(Arrow3D((x[0],y[0]),(x[1],y[1]),(x[2],y[2]),
          mutation_scale=11, lw=1, arrowstyle="-|>", color="g"))
        ax.text(x[0], x[1],x[2], iii, (1,0,0))
        ax.text(y[0], y[1],y[2], eee, (1,0,0))

# # # all face vectors.
#     for ix in faces_data['triangles_vertices']: #points3d:
#         vvv=ix
#         #print(vvv)
#         #ax.plot((vvv[0],), (vvv[1],), (vvv[2],), "o", color="#ff0000", ms=8, mew=0.5)
#         ax.add_artist(Arrow3D((0,vvv[0]),(0,vvv[1]),(0,vvv[2]),
#                               mutation_scale=1, lw=.5, arrowstyle="-", color="g"))

    ## Reset limits
#    ax.set_xlim(-1,1)
#    ax.set_ylim(-1,1)
#    ax.set_zlim(-1,1)
    ax.axis('off')
#    ax.view_init(elev=0, azim=60)
    ax.view_init(elev=0, azim=0)
    
    try:
        ax.set_aspect('equal')
    except NotImplementedError:
        pass
#        ax.set_box_aspect((1, 1, 1))
        
    ax.autoscale(True,axis='both')
    show()

        
if __name__ == "__main__":
    ##SC
    #faces_data = get_BZ(b1 = [1,0,0], b2 = [0,1,0], b3 = [0,0,1])
    ##BCC
    #faces_data = get_BZ(b1 = [1,1,0], b2 = [1,0,1], b3 = [0,1,1])
    ##FCC
    cell = [[0.5, 0.5, 1.0],[0.5, 1.0, 0.5],[1.0, 0.5, 0.5]]
#    cell=np.array([[8.885765876316732, 5.1301993206474545, 1.813799364234218], \
#          [-8.885765876316732, 5.1301993206474545, 1.813799364234218], \
#          [0.0, -10.260398641294914, 1.813799364234218]])/10.0
    qlat=[]
    vol=np.dot(cell[0],np.cross(cell[1],cell[2]))
    qlat.append(np.cross(cell[1],cell[2])/vol)
    qlat.append(np.cross(cell[2],cell[0])/vol)
    qlat.append(np.cross(cell[0],cell[1])/vol)
    cell=qlat
    plotws(cell[0],cell[1],cell[2])
