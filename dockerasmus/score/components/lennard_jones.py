# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

from .base import BaseComponent


class LennardJones(BaseComponent):
    """A scoring component modeling Van der Waals interactions.

    Reference:
        `Jones, J. E.
        "On the Determination of Molecular Fields. II. From the Equation of
        State of a Gas". Proceedings of the Royal Society of London A:
        Mathematical, Physical and Engineering Sciences 106, 463–477 (1924).
        <https://dx.doi.org/10.1098%2Frspa.1924.0082>`_
    """

    backends = ["theano", "mxnet", "tensorflow", "numpy"]

    def _setup_theano(self, theano):
        ### Potential well depth matrix from protein vectors
        v_pwd1, v_pwd2 = theano.tensor.dvectors('v_pwd1', 'v_pwd2')
        mx_well_depth = theano.tensor.sqrt(theano.tensor.outer(v_pwd1, v_pwd2))
        ### Radius matrix from protein vectors
        v_rad1, v_rad2 = theano.tensor.dvectors('v_rad1', 'v_rad2')
        mx_rad = theano.tensor.transpose(theano.tensor.repeat(
            theano.tensor.reshape(
                v_rad1, (1, v_rad1.shape[0])
            ), v_rad2.shape[0], axis=0,
        )) + v_rad2
        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')
        mx_distance_6 = mx_distance**6
        ### Van der Waals constants
        mx_rad_6 = mx_rad ** 6
        mx_B = 2 * mx_well_depth * mx_rad_6
        mx_A = 0.5 * mx_B * mx_rad_6
        ### Final function
        self._call = theano.function(
            [v_pwd1, v_pwd2, v_rad1, v_rad2, mx_distance],
            theano.tensor.sum(mx_A/(mx_distance_6**2) - mx_B/(mx_distance_6))
        )

    def _setup_numpy(self, numpy):
        from ...utils.matrices import compose
        def call(v_pwd1, v_pwd2, v_rad1, v_rad2, mx_distance):
            ### Potential well depth matrix from protein vectors
            mx_well_depth = numpy.sqrt(numpy.outer(v_pwd1, v_pwd2))
            ### Radius matrix from protein vectors (to the power of 6)
            mx_radius_6 = compose(lambda r1, r2: r1+r2, v_rad1, v_rad2)**6
            ### Van der Waals constants
            mx_B = 2 * mx_well_depth * mx_radius_6
            mx_A = 0.5 * mx_B * mx_radius_6
            ### Atomwise distance matrix
            mx_distance_6 = mx_distance**6
            ### Final function
            return numpy.sum(mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6))
        self._call = call

    def _setup_tensorflow(self, tf):
        tf_v_pwd1 = tf.placeholder(tf.float64)
        tf_v_pwd2 = tf.placeholder(tf.float64)
        tf_v_rad1 = tf.placeholder(tf.float64)
        tf_v_rad2 = tf.placeholder(tf.float64)

        mx_well_depth = tf.sqrt(tf.matmul(
            tf.expand_dims(tf_v_pwd1, 1),
            tf.expand_dims(tf_v_pwd2, 0),
        ))
        mx_radius_6 = (
            tf.transpose(tf.expand_dims(tf_v_rad1, 0)) + tf_v_rad2
        )**6
        mx_B = 2 * mx_well_depth * mx_radius_6
        mx_A = 0.5 * mx_B * mx_radius_6

        tf_mx_distance = tf.placeholder(tf.float64)
        mx_distance_6 = tf_mx_distance**6

        result = tf.reduce_sum(mx_A/(mx_distance_6**2) - mx_B/mx_distance_6)

        def call(v_pwd1, v_pwd2, v_rad1, v_rad2, mx_distance):
            with tf.Session() as sess:
                return sess.run(result, feed_dict={
                    tf_v_pwd1: v_pwd1,
                    tf_v_pwd2: v_pwd2,
                    tf_v_rad1: v_rad1,
                    tf_v_rad2: v_rad2,
                    tf_mx_distance: mx_distance,
                })

        self._call = call

        call([], [], [], [], [])


    def _setup_mxnet(self, mxnet):
        def call(v_pwd1, v_pwd2, v_rad1, v_rad2, mx_distance):
            ### Potential well depth matrix (sqrt of outer product)
            mx_well_depth = mxnet.ndarray.sqrt(mxnet.nd.dot(
                mxnet.nd.expand_dims(mxnet.nd.array(v_pwd1), 1),
                mxnet.nd.expand_dims(mxnet.nd.array(v_pwd2), 0),
            ))
            ### Radius matrix from protein vectors (to the power of 6)
            mx_radius_6 = (
                mxnet.nd.expand_dims(
                    mxnet.nd.array(v_rad1), 0,
                ).T + mxnet.nd.array(v_rad2)
            )**6
            ### Van der waals constants
            mx_B = 2 * mx_well_depth * mx_radius_6
            mx_A = 0.5 * mx_B * mx_radius_6
            ### Atomwise distance matrix (to the power of 6)
            mx_distance_6 = mxnet.nd.array(mx_distance)**6
            ### Final result
            return mxnet.nd.sum(
                mx_A/(mx_distance_6**2) - mx_B/mx_distance_6
            ).asscalar()
        self._call = call

    def __call__(self, potential_well_depth, vdw_radius, distance):
        return self._call(
            potential_well_depth[0],
            potential_well_depth[1],
            vdw_radius[0],
            vdw_radius[1],
            distance,
        )
