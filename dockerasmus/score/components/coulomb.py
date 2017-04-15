# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from .base import BaseComponent
from .. import requirements


class Coulomb(BaseComponent):
    """A scoring component modeling electrostatic forces.

    Reference:
        `Coulomb, C. A.
        "Premier mémoire sur l’électricité et le magnétisme". Histoire de
        l’Académie Royale des Sciences, 569-577 (1785).
        <https://books.google.com/books?id=by5EAAAAcAAJ&pg=PA569>`_
    """
    backends = ["theano", "tensorflow", "mxnet", "numpy"]


    def _setup_theano(self, theano):
        ### Dielectric constant
        diel = theano.tensor.dscalar('diel')
        ### Charge matrix from protein vectors
        v_q1, v_q2 = theano.tensor.dvectors('v_q1', 'v_q2')
        mx_q = theano.tensor.outer(v_q1, v_q2)
        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')
        ### Final function
        self._call = theano.function(
            [v_q1, v_q2, mx_distance, diel],
            theano.tensor.sum(mx_q/(diel*mx_distance))
        )

    def _setup_mxnet(self, mxnet):
        def call(v_q1, v_q2, mx_distance, diel):
            ### Convert input to mxnet NDArray
            v_q1 = mxnet.nd.array(v_q1)
            v_q2 = mxnet.nd.array(v_q2)
            mx_distance = mxnet.nd.array(mx_distance)
            ### Charge matrix from protein vectors
            mx_q = mxnet.nd.dot(
                v_q1.reshape((v_q1.shape[0],1)),
                v_q2.reshape((1, v_q2.shape[0]))
            )
            ### Atomwise distance matrix
            return mxnet.nd.sum(mx_q / (diel*mx_distance)).asscalar()
        self._call = call

    def _setup_tensorflow(self, tf):
        tf_v_q1 = tf.placeholder(tf.float64)
        tf_v_q2 = tf.placeholder(tf.float64)
        tf_mx_distance = tf.placeholder(tf.float64)
        tf_diel = tf.placeholder(tf.float64)

        tf_mx_q = tf.matmul(
            tf.expand_dims(tf_v_q1, -1),
            tf.expand_dims(tf_v_q2, 0),
        )
        result = tf.reduce_sum(tf_mx_q / (tf_diel * tf_mx_distance))

        def call(v_q1, v_q2, mx_distance, diel):
            with tf.Session() as sess:
                return sess.run(result, feed_dict={
                    tf_v_q1: v_q1,
                    tf_v_q2: v_q2,
                    tf_mx_distance: mx_distance,
                    tf_diel: diel,
                })

        self._call = call


    def _setup_numpy(self, numpy):
        def call(v_q1, v_q2, mx_distance, diel):
            mx_q = numpy.outer(v_q1, v_q2)
            return numpy.sum(mx_q/(diel*mx_distance))
        self._call = call

    def __call__(self, charge, distance, diel=65.0):
        return self._call(charge[0], charge[1], distance, diel)
