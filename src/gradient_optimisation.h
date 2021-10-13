/***************************************************************************
 *
 * Author: "Dari Kimanius"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef SOM_H
#define SOM_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <stack>
#include "src/parallel.h"
#include "src/filename.h"
#include <stdio.h>
#include <stdlib.h>

class SomGraph {

private:

	/**
	 * Class for graph nodes
	 */
	struct Node {
		bool active;
		float activity;
		float age;
	};

	std::vector<Node> _nodes;
	std::vector<float> _edges_activity;
	float neighbour_threshold;
	omp_lock_t mutex;

	/*
	 * Non-thread-safe add node.
	 */
	unsigned _add_node() {
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (!_nodes[i].active) { // If index not found
				_nodes[i].active = true;
				_nodes[i].activity = 0;
				_nodes[i].age = 0;
				return i;
			}
		throw std::runtime_error("failed to add node");
	}

	void _normalize_node_activity() {
		float activity_sum = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				activity_sum += _nodes[i].activity;

		if (activity_sum > 0)
			for (unsigned i = 0; i < _nodes.size(); i++)
				if (_nodes[i].active)
					_nodes[i].activity /= activity_sum;
	}

	void _normalize_edge_activity() {
		float activity_sum = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				for (unsigned j = i+1; j < _nodes.size(); j++)
					if (_nodes[j].active)
						activity_sum += _edges_activity[i * _nodes.size() + j];

		if (activity_sum > 0)
			for (unsigned i = 0; i < _edges_activity.size(); i++)
				_edges_activity[i] /= activity_sum;
	}

	/*
	 * Non-thread-safe remove node.
	 */
	void _remove_node(unsigned node) {
		if (!_nodes[node].active)
			throw std::runtime_error("node missing");

		_nodes[node].active = false;
		_nodes[node].activity = 0.;
		_nodes[node].age = 0.;

		for (unsigned i = 0; i < _nodes.size(); i++)
		{
			_edges_activity[node * _nodes.size() + i] = 0.;
			_edges_activity[i * _nodes.size() + node] = 0.;
		}

		_normalize_node_activity();
		_normalize_edge_activity();
	}

	unsigned _get_node_count() {
		unsigned count = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				count ++;
		return count;
	}

	void _init_mutex() {
		omp_init_lock(&mutex);
	}
public:

	SomGraph() {
		_init_mutex();
		neighbour_threshold = 0;
	}

	~SomGraph() {
		omp_destroy_lock(&mutex);
	}

	/**
	 * Set max number of nodes
	 */
	void set_max_node_count(unsigned count) {
		Node n;
		n.active = false;
		n.activity = 0;
		n.age = 0.;

		_nodes.resize(count, n);
		_edges_activity.resize(count * count, 0);
	}

	/**
	 * Clone nodes settings
	 */
	void clone_nodes(const SomGraph &som) {
		_nodes = som._nodes;
		_edges_activity.resize(som._edges_activity.size(), 0);
	}

	void reset_activities() {
		for (unsigned i = 0; i < _nodes.size() + 1; i++)
			_nodes[i].activity = 0.;
		for (unsigned i = 0; i < _edges_activity.size() + 1; i++)
			_edges_activity[i] = 0.;
	}

	/**
	 * Add an unconnected node to the graph.
	 */
	unsigned add_node() {
		Lock ml(&mutex);
		return _add_node();
	}

	/**
	 * Add a node to the graph close to given node.
	 */
	unsigned add_node(unsigned node, float age_factor=0.1) {
		Lock ml(&mutex);
		unsigned n = _add_node();
		_edges_activity[node * _nodes.size() + n] = neighbour_threshold;
		_edges_activity[n * _nodes.size() + node] = neighbour_threshold;
		_nodes[n].activity = _nodes[node].activity;
		_nodes[n].age = _nodes[node].age * age_factor;
		_normalize_node_activity();
		_normalize_edge_activity();
		return n;
	}

	/**
	 * Remove node.
	 */
	void remove_node(unsigned node) {
		Lock ml(&mutex);
		_remove_node(node);
	}

	/**
	 * Getters and setters.
	 */

	float get_edge_activity(unsigned node1, unsigned node2) {
		if (node1 == node2)
			throw std::runtime_error("edge to node itself");
		Lock ml(&mutex);
		return _edges_activity[node1 * _nodes.size() + node2];
	}

	void set_edge_activity(unsigned node1, unsigned node2, float activity=1.) {
		if (node1 == node2)
			throw std::runtime_error("edge to node itself");
		Lock ml(&mutex);
		_edges_activity[node1 * _nodes.size() + node2] = activity;
		_edges_activity[node2 * _nodes.size() + node1] = activity;
	}

	void add_edge_activity(unsigned node1, unsigned node2, float activity=1.) {
		if (node1 == node2)
			throw std::runtime_error("edge to node itself");
		Lock ml(&mutex);
		_edges_activity[node1 * _nodes.size() + node2] += activity;
		_edges_activity[node2 * _nodes.size() + node1] += activity;
	}

	float get_node_age(unsigned node) {
		Lock ml(&mutex);
		return _nodes[node].age;
	}

	void set_node_age(unsigned node, float age=0.) {
		Lock ml(&mutex);
		_nodes[node].age = age;
	}

	void add_node_age(unsigned node, float age=1.) {
		Lock ml(&mutex);
		_nodes[node].age += age;
	}

	float get_node_activity(unsigned node) {
		Lock ml(&mutex);
		return _nodes[node].activity;
	}

	void set_node_activity(unsigned node, float activity=0.) {
		Lock ml(&mutex);
		_nodes[node].activity = activity;
	}

	void add_node_activity(unsigned node, float activity=0.) {
		Lock ml(&mutex);
		_nodes[node].activity += activity;
	}

	unsigned get_node_count() {
		Lock ml(&mutex);
		return _get_node_count();
	}

	void normalize_activity() {
		_normalize_node_activity();
		_normalize_edge_activity();
	}

	std::vector<unsigned> get_all_nodes() {
		Lock ml(&mutex);
		std::vector<unsigned> nodes;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				nodes.push_back(i);
		return nodes;
	}

	void update_node_activities(SomGraph update, float mu) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _nodes.size(); i++)
			_nodes[i].activity = _nodes[i].activity * mu + update._nodes[i].activity * (1-mu);
		_normalize_node_activity();
	}

	void update_edge_activities(SomGraph update, float mu) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges_activity.size(); i++)
			_edges_activity[i] = _edges_activity[i] * mu + update._edges_activity[i] * (1-mu);
		_normalize_edge_activity();
	}

	std::vector< std::pair<unsigned, float> > get_neighbours(unsigned node) {
		Lock ml(&mutex);
		std::vector< std::pair<unsigned, float> > out;
		float weight_sum = 0;
		for (unsigned i = 0; i < _nodes.size(); i++) {
			if (i == node)
				continue;
			float w = _edges_activity[_nodes.size() * node + i];
			if (w > neighbour_threshold) {
				out.push_back(std::pair<unsigned, float>(i, w));
				weight_sum += w;
			}
		}

		for (unsigned i = 0; i < out.size(); i++)
			out[i].second /= weight_sum;

		return out;
	}

	/*
	 * Set the averactivity number of neighbours per node.
	 */
	void set_connectivity(float connectivity) {
		if (connectivity < 0)
			throw std::runtime_error("bad connectivity value");

		Lock ml(&mutex);
		_normalize_edge_activity();

		unsigned N = _get_node_count();
		std::vector<float> edges(( N * (N-1) ) / 2, 0);

		int n = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				for (unsigned j = i+1; j < _nodes.size(); j++)
					if (_nodes[j].active) {
						edges[n] = _edges_activity[i * _nodes.size() + j];
						n ++;
					}
		std::sort (edges.begin(), edges.end());
		int nr_edges = N * connectivity / 2;  // Divide by 2, each edge is two connections
		if (nr_edges >= edges.size())
			neighbour_threshold = 0;
		else
			neighbour_threshold = edges[edges.size() - nr_edges];
	}

	void print_to_file(FileName &fn) {
		Lock ml(&mutex);
		FILE * fp;
		fp = fopen (fn.c_str(), "w+");

		fprintf(fp, "\n_edge_activity_threshold %4.4f\n", neighbour_threshold);

		fprintf(fp, "\n_nodes [index age activity]\n");

		std::vector<int> active_nodes;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active) {
				fprintf(fp, "%3d %4.4f %4.4f\n", i, _nodes[i].age, _nodes[i].activity);
				active_nodes.push_back(i);
			}

		fprintf(fp, "\n_edge_activities\n");

		for(unsigned i = 0; i < active_nodes.size(); i++) {
			for (unsigned j = 0; j < active_nodes.size(); j++)
				fprintf(fp, "%4.4f ", _edges_activity[active_nodes[i] * _nodes.size() + active_nodes[j]]);
			fprintf(fp, "\n");
		}

		fclose(fp);
	}

	template<typename T>
	static bool _pair_sort_ascend_fn (std::pair<unsigned, T> i, std::pair<unsigned, T> j)
	{ return (i.second < j.second); }

	template<typename T>
	static bool _pair_sort_descend_fn (std::pair<unsigned, T> i, std::pair<unsigned, T> j)
	{ return (i.second > j.second); }

	template<typename T>
	static std::vector<unsigned> arg_sort(std::vector<T> values, bool ascend=true) {
		std::vector< std::pair<unsigned, T> > pair_list(values.size());
		for (unsigned i = 0; i < values.size(); i ++) {
			pair_list[i].first = i;
			pair_list[i].second = values[i];
		}
		if (ascend)
			std::sort (pair_list.begin(), pair_list.end(), _pair_sort_ascend_fn<T>);
		else
			std::sort (pair_list.begin(), pair_list.end(), _pair_sort_descend_fn<T>);

		std::vector<unsigned> list(values.size());
		for (unsigned i = 0; i < values.size(); i ++) {
			list[i] = pair_list[i].first;
		}
		return list;

	}

	// Generate random number on a normal distribution, with std=1 and mean=0
	template<typename T>
	static T random_normal()
	{
		T v1 = ((T) (rand()) + 1.) / ((T) (RAND_MAX) + 1.);
		T v2 = ((T) (rand()) + 1.) / ((T) (RAND_MAX) + 1.);
		return cos(2 * 3.14 * v2) * sqrt(-2. * log(v1));
	}

	template<typename T>
	static T mean(MultidimArray<T> &array)
	{
		return array.sum() / array.nzyxdim;
	}

	template<typename T>
	static T variance(MultidimArray<T> &array)
	{
		T diff, std(0), avg(mean(array));
		for (int i = 0; i < array.nzyxdim; i ++)
		{
			diff = array.data[i] - avg;
			std += diff * diff;
		}
		std /= (T) array.nzyxdim;
		return std;
	}

	template<typename T>
	static T std(MultidimArray<T> &array)
	{
		return sqrt(variance(array));
	}

	template<typename T>
	static void make_blobs_2d(MultidimArray<T> &box, MultidimArray<T> &amp_box, unsigned nr_blobs, T diameter, bool helical)
	{
		std::vector<T> blobs_x(nr_blobs), blobs_y(nr_blobs), blobs_amp(nr_blobs, 0);
		box.resize(amp_box);
		box.initZeros();

		for (int i = 0; i < nr_blobs; i ++)
		{
			// Generate random coordinates
			if (helical)
				blobs_x[i] = (T) (rand() % box.xdim);
			else
				blobs_x[i] = random_normal<T>() * diameter / 7. + box.xdim / 2.;
			blobs_y[i] = random_normal<T>() * diameter / 7. + box.ydim / 2.;
		}

		for (int i = 0; i < nr_blobs; i ++)
			blobs_amp[i] = fabs(amp_box.data[(long) blobs_y[i] * amp_box.xdim + (long) blobs_x[i]]);

		T sigma_inv = 10. / diameter;
		for (int y = 0; y < box.ydim; y ++)
			for (int x = 0; x < box.xdim; x ++)
				for (int i = 0; i < nr_blobs; i ++)
				{
					T xp = (x-blobs_x[i]) * sigma_inv;
					T yp = (y-blobs_y[i]) * sigma_inv;
					box.data[y * box.xdim + x] += exp(-xp*xp -yp*yp) * blobs_amp[i];
				}
	}

	template<typename T>
	static void make_blobs_3d(MultidimArray<T> &box, MultidimArray<T> &amp_box, unsigned nr_blobs, T diameter, bool helical)
	{
		std::vector<T> blobs_x(nr_blobs), blobs_y(nr_blobs), blobs_z(nr_blobs), blobs_amp(nr_blobs);
		box.resize(amp_box);
		box.initZeros();

		for (int i = 0; i < nr_blobs; i ++)
		{
			// Generate random coordinates
			blobs_x[i] = random_normal<T>() * diameter / 6. + box.xdim / 2.;
			blobs_y[i] = random_normal<T>() * diameter / 6. + box.ydim / 2.;

			if (helical)
				blobs_z[i] = ((T) rand() / (T) RAND_MAX) * box.zdim;
			else
				blobs_z[i] = random_normal<T>() * diameter / 6. + box.zdim / 2.;
		}

		for (int i = 0; i < nr_blobs; i ++)
			blobs_amp[i] = fabs(amp_box.data[
				(long) blobs_z[i] * amp_box.yxdim +
				(long) blobs_y[i] * amp_box.xdim +
				(long) blobs_x[i]
			]);

		T span = diameter/3.;
		T sigma_inv = 10./diameter;
		for (int i = 0; i < nr_blobs; i ++)
			for (int z = XMIPP_MAX(0, blobs_z[i]-span);
				z < XMIPP_MIN(box.zdim, blobs_z[i]+span); z ++)
				for (int y = XMIPP_MAX(0, blobs_y[i]-span);
					y < XMIPP_MIN(box.ydim, blobs_y[i]+span); y ++)
					for (int x = XMIPP_MAX(0, blobs_x[i]-span);
						x < XMIPP_MIN(box.xdim, blobs_x[i]+span); x ++)
					{
						T xp = (x-blobs_x[i]) * sigma_inv;
						T yp = (y-blobs_y[i]) * sigma_inv;
						T zp = (z-blobs_z[i]) * sigma_inv;
						box.data[z * amp_box.yxdim + y * box.xdim + x] +=
								exp(-xp*xp -yp*yp -zp*zp) * blobs_amp[i];
					}
	}

};

#endif
