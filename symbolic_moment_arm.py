# \brief Calculates a symbolic expression of the muscle moment arm of
#  an OpenSim .osim model. The moment arm is sampled and approximated
#  by a multivariate polynomial, so that higher order derivatives can
#  be computed. This implementation works with OpenSim v4.0 API.
#
# Dependencies: opensim, matplotlib, numpy, sympy, multipolyfit
#
# @author Dimitar Stanev (jimstanev@gmail.com)
import csv
import pickle
import opensim
import collections
import numpy as np
import sympy as sp
import operator # used in sorted
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # for projection='3d'
from multipolyfit import multipolyfit, mk_sympy_function


###############################################################################
# parameters

model_file = 'Hamner2010_patella_scaled.osim'
muscle_coordinates_file = 'muscle_coordinates.csv'
# when computed once results are stored into files and loaded with (pickle)
compute = False

###############################################################################
# utilities

def cartesian(arrays, out=None):
    """Generate a cartesian product of input arrays.

    Parameters
    ----------

    arrays: list of array-like
        1-D arrays to form the cartesian product of.

    out: ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out: ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])
    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m, 1:] = out[0:m, 1:]

    return out


def construct_coordinate_grid(model, coordinates, N=5):
    """Given n coordinates get the coordinate range and generate a coordinate grid
    of combinations using cartesian product.

    Parameters
    ----------

    model: opensim.Model

    coordinates: list of string

    N: int (default=5)
        the number of points per coordinate

    Returns
    -------

    sampling_grid: np.array
        all combination of coordinates

    """
    sampling_grid = []
    for coordinate in coordinates:
        min_range = model.getCoordinateSet().get(coordinate).getRangeMin()
        max_range = model.getCoordinateSet().get(coordinate).getRangeMax()
        sampling_grid.append(np.linspace(min_range, max_range, N,
                                         endpoint=True))

    return cartesian(sampling_grid)


def visualize_moment_arm(moment_arm_coordinate, muscle, coordinates,
                         sampling_dict, model_coordinates, model_muscles, R):
    """Visualize moment arm as 2D or 3D plot.

    Parameters
    ----------

    moment_arm_coordinate: string
        which moment arm (coordinate)

    muscle: string
        which muscle

    coordinates: list of strings
        which coordinates affect the moment arm variable (one or two only)

    sampling_dict: dictionary
        calculated from calculate_moment_arm_symbolically

    model_coordinates: dictionary
        coordinate names and their corresponding indices in the model

    model_muscles: dictionary
        muscle names and their corresponding indices in the model
    """
    if isinstance(coordinates, str):
        coordinates = sampling_dict[muscle]['coordinates']
        sampling_grid = sampling_dict[muscle]['sampling_grid']
        moment_arm = sampling_dict[muscle]['moment_arm']
        idx = coordinates.index(moment_arm_coordinate)
        poly = R[model_muscles[muscle],
                 model_coordinates[moment_arm_coordinate]]
        moment_arm_poly = np.array([poly.subs(dict(zip(poly.free_symbols, x)))
                                    for x in sampling_grid],
                                   np.float)

        fig = plt.figure()
        ax = fig.gca()
        ax.plot(sampling_grid[:, idx], moment_arm[:, idx], 'rx',
                label='sampled')
        ax.plot(sampling_grid[:, idx], moment_arm_poly, 'b-',
                label='fitted')
        ax.set_xlabel(coordinates[idx])
        ax.set_ylabel(muscle + '_' + moment_arm_coordinate)
        ax.set_title('2D Moment Arm')
        ax.legend()
        fig.tight_layout()
        fig.savefig(muscle + '_' + moment_arm_coordinate + '.pdf',
                    format='pdf', dpi=300)
        fig.savefig(muscle + '_' + moment_arm_coordinate + '.png',
                    format='png', dpi=300)
        plt.show()
    elif isinstance(coordinates, list) and len(coordinates) == 2:
        coordinates = sampling_dict[muscle]['coordinates']
        sampling_grid = sampling_dict[muscle]['sampling_grid']
        moment_arm = sampling_dict[muscle]['moment_arm']
        idx = coordinates.index(moment_arm_coordinate)
        poly = R[model_muscles[muscle],
                 model_coordinates[moment_arm_coordinate]]
         # poly.free_symbols is not used because it may not preserve order
        poly_symbols = [sp.Symbol(x) for x in coordinates]
        moment_arm_poly = np.array([
            poly.subs(dict(zip(poly_symbols, x))) for x in sampling_grid
        ], np.float)

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(sampling_grid[:, 0], sampling_grid[:, 1],
                   moment_arm[:, idx],
                   label='sampled', color='r')
        surf = ax.plot_trisurf(sampling_grid[:, 0], sampling_grid[:, 1],
                               moment_arm_poly, label='fitted',
                               facecolor='b', edgecolor='k',
                               linewidth=0.1, antialiased=True)
        surf._facecolors2d = surf._facecolors3d
        surf._edgecolors2d = surf._edgecolors3d
        ax.set_xlabel(coordinates[0])
        ax.set_ylabel(coordinates[1])
        ax.set_zlabel(muscle + '_' + moment_arm_coordinate)
        ax.set_title('3D Moment Arm')
        ax.legend()
        fig.tight_layout()
        fig.savefig(muscle + '_' + moment_arm_coordinate + '.pdf',
                    format='pdf', dpi=300)
        fig.savefig(muscle + '_' + moment_arm_coordinate + '.png',
                    format='png', dpi=300)
        plt.show()
    else:
        return


def calculate_moment_arm_symbolically(model_file, muscle_coordinates_file):
    """Calculate the muscle moment arm matrix symbolically for a particular
     OpenSim model.

    Parameters
    ----------

    model_file: string
        OpenSim .osim

    muscle_coordinates_file: string
        a csv file (';' delimiter) containing the name of the muscle
        and the spanning coordinates in each row
    """
    # parse csv
    muscle_coordinates = {}
    with open(muscle_coordinates_file) as csv_file:
        reader = csv.reader(csv_file, delimiter=';')
        for row in reader:
            muscle_coordinates[row[0]] = row[1:]

    # load opensim model
    model = opensim.Model(model_file)
    state = model.initSystem()

    model_coordinates = {}
    for i in range(0, model.getNumCoordinates()):
        model_coordinates[model.getCoordinateSet().get(i).getName()] = i

    model_muscles = {}
    for i in range(0, model.getNumControls()):
        model_muscles[model.getMuscles().get(i).getName()] = i

    # calculate moment arm matrix (R) symbolically
    R = []
    sampling_dict = {}
    resolution = {1: 15, 2: 10, 3: 8, 4: 5, 5: 5}
    for muscle, k in sorted(model_muscles.items(), key=operator.itemgetter(1)):
        # get initial state each time
        state = model.initSystem()
        coordinates = muscle_coordinates[muscle]
        N = resolution[len(coordinates)]
        print(N, muscle, coordinates)
        # calculate moment arms for this muscle and spanning coordinates
        sampling_grid = construct_coordinate_grid(model, coordinates, N)
        moment_arm = []
        for q in sampling_grid:
            for i, coordinate in enumerate(coordinates):
                model.updCoordinateSet().get(coordinate).setValue(state, q[i])

            model.realizePosition(state)
            tmp = []
            for coordinate in coordinates:
                coord = model.getCoordinateSet().get(coordinate)
                tmp.append(model.getMuscles()
                           .get(muscle).computeMomentArm(state, coord))

            moment_arm.append(tmp)

        moment_arm = np.array(moment_arm)
        sampling_dict[muscle] = {'coordinates': coordinates,
                                 'sampling_grid': sampling_grid,
                                 'moment_arm': moment_arm}

        # polynomial regression
        degree = 5
        muscle_moment_row = [0] * len(model_coordinates)
        for i, coordinate in enumerate(coordinates):
            coeffs, powers = multipolyfit(sampling_grid, moment_arm[:, i],
                                          degree, powers_out=True)
            polynomial = mk_sympy_function(coeffs, powers)
            polynomial = polynomial.subs(dict(zip(polynomial.free_symbols,
                                                  [sp.Symbol(x) for x in
                                                   coordinates])))
            muscle_moment_row[model_coordinates[coordinate]] = polynomial

        R.append(muscle_moment_row)

    # export data to file because the process is time consuming
    R = sp.Matrix(R)
    pickle.dump(R, file('R.dat', 'w'))
    pickle.dump(sampling_dict, file('sampling_dict.dat', 'w'))
    pickle.dump(model_muscles, file('model_muscles.dat', 'w'))
    pickle.dump(model_coordinates, file('model_coordinates.dat', 'w'))


def find_intermediate_joints(origin_body, insertion_body, model_tree, joints):
    """Finds the intermediate joints between two bodies.

    Parameters
    ----------

    origin_body: string
        first body in the model tree

    insertion_body: string
        last body in the branch

    model_tree: list of dictionary relations {parent, joint, child}

    joints: list of strings
        intermediate joints
    """
    if origin_body == insertion_body:
        return True

    children = filter(lambda x: x['parent'] == origin_body, model_tree)
    for child in children:
        found = find_intermediate_joints(child['child'], insertion_body,
                                         model_tree, joints)
        if found:
            joints.append(child['joint'])
            return True

    return False


def calculate_spanning_muscle_coordinates(model_file, muscle_coordinates_file):
    """Calculates the coordinates that are spanned by each muscle. Useful for
    reducing the required computation of the muscle moment arm matrix.

    Parameters
    ----------

    model_file: string
        OpenSim .osim

    muscle_coordinates_file: string
        where to store the results as a csv file (';' delimiter) containing
        the name of the muscle and the spanning coordinates in each row

    """
    model = opensim.Model(model_file)
    model.initSystem()

    # construct model tree (parent body - joint - child body)
    model_tree = []
    for joint in model.getJointSet():
        model_tree.append({'parent': joint.getParentFrame().getName()
                           .replace('_offset', ''), # v4.0 convention
                           'joint': joint.getName(),
                           'child': joint.getChildFrame().getName()})

    # get the coordinates that are spanned by the muscles
    muscle_coordinates = {}
    for muscle in model.getMuscles():
        path = muscle.getGeometryPath().getPathPointSet()
        muscle_bodies = []
        for point in path:
            muscle_bodies.append(point.getBodyName())

        # get unique bodies (e.g. as in set()) and preserve order of insertion
        muscle_bodies = collections.OrderedDict.fromkeys(muscle_bodies).keys()

        # find intermediate joints
        assert(len(muscle_bodies) > 1)
        joints = []
        find_intermediate_joints(muscle_bodies[0], muscle_bodies[-1],
                                 model_tree, joints)

        # find spanning coordinates
        muscle_coordinates[muscle.getName()] = []
        for joint in joints:
            joint = model.getJointSet().get(joint)
            for i in range(0, joint.numCoordinates()):
                muscle_coordinates[muscle.getName()].append(
                    joint.get_coordinates(i).getName())

    # write results to file
    with open(muscle_coordinates_file, 'w') as csv_file:
        for key, values in muscle_coordinates.items():
            csv_file.write(key)
            for value in values:
                csv_file.write(';' + value)

            csv_file.write('\n')


###############################################################################
# main

if compute:
    calculate_spanning_muscle_coordinates(model_file, muscle_coordinates_file)
    calculate_moment_arm_symbolically(model_file, muscle_coordinates_file)

# load data
R = pickle.load(file('R.dat', 'r'))
sampling_dict = pickle.load(file('sampling_dict.dat', 'r'))
model_coordinates = pickle.load(file('model_coordinates.dat', 'r'))
model_muscles = pickle.load(file('model_muscles.dat', 'r'))

# visualize 3D moment arm R(q1, q2)
muscle = 'tib_ant_r'
coordinates = sampling_dict[muscle]['coordinates']
visualize_moment_arm(coordinates[1], muscle, coordinates, sampling_dict,
                     model_coordinates, model_muscles, R)

# visualize 2D moment arm R(q1)
muscle = 'vas_int_r'
coordinates = sampling_dict[muscle]['coordinates']
visualize_moment_arm(coordinates[0], muscle, coordinates[0], sampling_dict,
                     model_coordinates, model_muscles, R)
