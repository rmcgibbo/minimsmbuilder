import os, sys
import numpy as np
import random

# Serializer is a convenience wrapper around the python HDF5 library
# for reading in data from disk (and saving results)
from minimsmbuilder.Serializer import Serializer

# Project subclasses serializer, and is used to read in a file called
# "ProjectInfo".h5 that contains data about things like the path to each of
# the trajectory objects and such
from minimsmbuilder.Project import Project

# the Trajectory object is basically a dictionary. The main thing it contains
# is "XYZList", and N x M x 3 matrix, where N is the number of protein conformations
# in the trajectory, M is the number of atoms in the protein, and 3 for each atoms'
# X, Y and Z coordinates
from minimsmbuilder.Trajectory import Trajectory

# The TheoData object handles preparing the xyz coordinates for the RMSD calculation
# in C. This preparation involves making sure everything is byte-aligned and padded with
# the appropriate zeros, translating each conformation so that its center of mass
# is at the origin, etc.
from minimsmbuilder.TheoData import TheoData

# This is the C extension code that operated on the TheoData objects
from minimsmbuilder import rmsdcalc


def one_to_all(Theo1, Theo2, index1):
    '''Calculate a vector of distances from the ith frame of Theo1 to all the frames
    of Theo2.
    
    Arguments:
    Theo1 - TheoData object
    Theo2 - TheoData object
    index1 - integer index into TheoData1
    
    Returns:
    Vector of floats of length equal to the length of Theo2.
    '''
    return rmsdcalc.getMultipleRMSDs_aligned_T_g(
        Theo1.NumAtoms, Theo1.NumAtomsWithPadding,
        Theo1.NumAtomsWithPadding, Theo2.XYZData,
        Theo1.XYZData[index1], Theo2.G,
        Theo1.G[index1])


def assign(ptraj, generator_indices):
    '''Given a TheoData (ptraj) and some indices to be the centers,
    compute the assignments and distances for each of the data points
    '''
    assignments = np.zeros(len(ptraj), dtype='int')
    distances = np.inf * np.ones(len(ptraj), dtype='float32')
    for m in generator_indices:
        d = one_to_all(ptraj, ptraj, m)
        closer = np.where(d < distances)[0]
        distances[closer] = d[closer]
        assignments[closer] = m
    return assignments, distances


def kcenters(theo, k):
    '''
    Cluster the conformations in the prepared trajectory (Theo) into k clusters based on the RMSD distance
    metric and the kcenters clustering algorithm.
    
    This is a cheap (and pretty bad) clustering algorithm
    '''
    
    number_of_conformations = len(theo.XYZData)
    
    new_center = 0 # we start by selecting (arbitrarily) the zeroth frame to be
    # initial center
    center_indices = [] #indices of the identified kcenters
    # distance_list holds the distance of every conformation to its closest center
    distance_list = np.inf * np.ones(number_of_conformations)
    # assignments holds the index of the center that each conformation is closest to
    assignments = -1 * np.ones(number_of_conformations, dtype='int')
    
    for i in xrange(k):
        print 'Finding Generator %d' % i
        center_indices.append(new_center)
        distance_to_new_center = one_to_all(theo, theo, new_center)
        updated_indices = np.where(distance_to_new_center < distance_list)[0]
        distance_list[updated_indices] = distance_to_new_center[updated_indices]
        assignments[updated_indices] = new_center
        
        # pick the data point that's currently worst assigned as our new center
        new_center = np.argmax(distance_list)
    
    return np.array(center_indices), assignments, distance_list


def clarans(ptraj, k, num_local_minima, max_neighbors, local_swap=True, initial_medoids='kcenters', initial_assignments=None, initial_distance=None, verbose=True):
    '''
    Cluster the conformations in the prepared trajectory (Theo) into k clusters based on the RMSD distance
    metric and the CLARANS clustering algorithm.
    
    http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=01033770
    
    Accuracy and cost are both higher as you increase the num_local_minima and max_neighbors parametors.
    The other parameters should probably be left with their default values.
    '''
    
    
    num_frames = len(ptraj)

    if initial_medoids == 'kcenters':
        initial_medoids, initial_assignments, initial_distance = kcenters(ptraj, k)
    elif initial_medoids == 'random':
        initial_medoids = np.random.permutation(np.arange(num_frames))[0:k]
        initial_assignments, initial_distance = assign(ptraj, initial_medoids)
    else:
        if not isinstance(initial_medods, np.ndarray):
            raise ValueError('Initial medoids should be a numpy array')
        if initial_assignments is None or initial_distance is None:
            initial_assignments, initial_distance = assign(ptraj, initial_medoids)

    if not len(initial_assignments) == num_frames:
        raise ValueError('Initial assignments is not the same length as ptraj')
    if not len(initial_distance) == num_frames:
        raise ValueError('Initial distance is not the same length as ptraj')
    if not k == len(initial_medoids):
        raise ValueError('Initial medoids not the same length as k')

    if verbose:
        printer = sys.stdout
    else:
        printer = open('/dev/null', 'w')


    initial_pmedoids = ptraj[initial_medoids]
    initial_cost = np.sum(initial_distance)
    min_cost = initial_cost

    # these iterations could be parallelized
    for i in xrange(num_local_minima):
        print >> printer, '%s of %s local minima' % (i, num_local_minima)

        # the cannonical clarans approach is to initialize the medoids that you
        # start from randomly, but instead we use the kcenters medoids.

        medoids = initial_medoids
        pmedoids = initial_pmedoids
        assignments = initial_assignments
        distance_to_current = initial_distance
        current_cost = initial_cost


        #loop over neighbors
        j = 0
        while j < max_neighbors:
            medoid_i = np.random.randint(k)
            old_medoid = medoids[medoid_i]

            if local_swap is False:
                trial_medoid = np.random.randint(num_frames)
            else:
                trial_medoid = random.choice(np.where(assignments == medoids[medoid_i])[0])

            new_medoids = medoids.copy()
            new_medoids[medoid_i] = trial_medoid
            pmedoids = ptraj[new_medoids]

            new_distances = distance_to_current.copy()
            new_assignments = assignments.copy()

            print >> printer, '  swapping %s for %s...' % (old_medoid, trial_medoid),

            distance_to_trial = one_to_all(ptraj, ptraj, trial_medoid)
            assigned_to_trial = np.where(distance_to_trial < distance_to_current)[0]
            new_assignments[assigned_to_trial] = trial_medoid
            new_distances[assigned_to_trial] = distance_to_trial[assigned_to_trial]

            ambiguous = np.where((new_assignments == old_medoid) & \
                                 (distance_to_trial >= distance_to_current))[0]
            for l in ambiguous:
                d = one_to_all(ptraj, pmedoids, l)
                argmin =  np.argmin(d)
                new_assignments[l] = new_medoids[argmin]
                new_distances[l] = d[argmin]

            new_cost = np.sum(new_distances)
            if new_cost < current_cost:
                print >> printer, 'Accept'
                medoids = new_medoids
                assignments = new_assignments
                distance_to_current = new_distances
                current_cost = new_cost

                j = 0
            else:
                j += 1
                print >> printer, 'Reject'

        if current_cost < min_cost:
            min_cost = current_cost
            optimal_medoids = medoids.copy()
            optimal_assignments = assignments.copy()
            optimal_distances = distance_to_current.copy()


    return optimal_medoids, optimal_assignments, optimal_distances


def main():
    data_dir = os.path.join(os.path.dirname(__file__), 'data/')
    project_info_path = os.path.join(data_dir, 'ProjectInfo.h5')

    project = Project.LoadFromHDF(project_info_path)
    trajectory_xyzlist = project.GetAllConformations()

    # If you want to simulate a longer trajectory (these trajectories)
    # are small compared to our real datasets, you can just stack this
    # one on top of itself:
    # multiplier = 5
    # trajectory_xyzlist = np.vstack([trajectory_xyzlist] * multiplier)

    # The TheoData object prepares the xyz coordinates for the C code. This
    # preparation involves byte-aligning some of the arrays, translating 
    # the xyz coordinates of each conformation so that its center of mass is at
    # the origin, and precalculating some value called G (some kind of magnitude)
    # for each of the conformations, which is required by the RMSD algorithm.
    theo = TheoData(trajectory_xyzlist)

    num_clusters = 100

    # run clarans clustering
    # The execution time is mostly set by num_local_minima and max_neighbors
    # (lower numbers go faster). The execution time should be roughly
    # to num_local_minima, and scales worse with max_neighbors
    #
    # Also, as you scale to bigger num_local_minima and max_neighbors, the percentage
    # of the execution time spent in the C extension code goes up.
    centers, assignments, distances = clarans(theo, num_clusters,
                                              num_local_minima=5, max_neighbors=5)


    print '\nFound %d centers with kcenters' % num_clusters
    print centers
    
if __name__ == '__main__':
    main()

