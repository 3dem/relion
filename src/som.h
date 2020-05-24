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
#include <unordered_map>
#include <unordered_set>
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
		bool active = false;
		float age = 0;
		float error = 0.;
	};

	/**
	 * Class for graph edges
	 */
	class Edge {
	public:
		float age = 0.;
		unsigned n1, n2;
		Edge (unsigned node1, unsigned node2):
				n1(node1), n2(node2)
		{};
	};

	std::vector<Node> _nodes;
	std::vector<Edge> _edges;

	pthread_mutex_t mutex;

	/*
	 * Non-thread-safe remove node.
	 */
	void _remove_node(unsigned node) {
		if (!_nodes[node].active)
			throw std::runtime_error("node missing");

		_nodes[node].active = false;

		for (unsigned i = 0; i < _edges.size(); i++)
			if (_edges[i].n1 == node || _edges[i].n2 == node)
				_edges.erase(_edges.begin()+i);

	}

	/*
	 * Non-thread-safe get neighbours of node.
	 */
	std::vector<unsigned> _get_neighbours(unsigned node) {
		std::vector<unsigned> neighbours;
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node)
				neighbours.push_back(_edges[i].n2);
			if (_edges[i].n2 == node)
				neighbours.push_back(_edges[i].n1);
		}
		return neighbours;
	}

public:

	SomGraph()
	{
		int mutex_error = pthread_mutex_init(&mutex, NULL);

		if (mutex_error != 0)
			throw std::runtime_error("failed to initialize mutex");
	}

	/**
	 * Max number of nodes
	 */
	void set_max_node_count(unsigned count) {
		_nodes.resize(count);
	}

	/**
	 * Add an edge-less node to the graph.
	 */
	unsigned add_node() {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _nodes.size() + 1; i++)
			if (!_nodes[i].active) { // If index not found
				_nodes[i].active = true;
				_nodes[i].age = 0;
				_nodes[i].error = 0;
				return i;
			}
		throw std::runtime_error("failed to add node");
	}

	/**
	 * Add a edge between node1 and node2.
	 * Do nothing if edge already exists.
	 */
	void add_edge(unsigned node1, unsigned node2) {
		if (node1 == node2)
			throw std::runtime_error("cannot add edge to the same node itself");

		if (!_nodes[node1].active || !_nodes[node2].active)
			throw std::runtime_error("node(s) missing");

		Lock ml(&mutex);

		// Does edge already exist
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node1 && _edges[i].n2 == node2 ||
			    _edges[i].n1 == node2 && _edges[i].n2 == node1)
				return;
		}

		_edges.emplace_back(node1, node2);
	}

	/**
	 * Remove node.
	 */
	void remove_node(unsigned node) {
		Lock ml(&mutex);
		_remove_node(node);
	}

	/**
	 * Remove edge.
	 */
	void remove_edge(unsigned node1, unsigned node2) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node1 && _edges[i].n2 == node2 ||
			    _edges[i].n1 == node2 && _edges[i].n2 == node1) {
				_edges.erase(_edges.begin()+i);
				return;
			}
		}
		throw std::runtime_error("edge not found");
	}

	/**
	 * Get neighbours of given node.
	 */
	std::vector<unsigned> get_neighbours(unsigned node) {
		Lock ml(&mutex);
		return _get_neighbours(node);
	}

	/**
	 * Get average edge age.
	 */
	float get_avg_age() {
		Lock ml(&mutex);
		float avg = 0;
		for (unsigned i = 0; i < _edges.size(); i++)
			avg += _edges[i].age;
		return avg / (float) _edges.size();
	}

	/**
	 * Get average node error.
	 */
	float get_avg_error() {
		Lock ml(&mutex);
		float avg = 0;
		int count = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active) {
				avg += _nodes[i].error;
				count ++;
			}
		return avg / count;
	}

	/**
	 * Get neighbours of given node and their corresponding age.
	 */
	std::vector<std::pair<unsigned, float>> get_neighbours_age(unsigned node) {
		Lock ml(&mutex);
		std::vector<std::pair<unsigned, float>> out;
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node)
				out.emplace_back(_edges[i].n2, _edges[i].age);
			if (_edges[i].n2 == node)
				out.emplace_back(_edges[i].n1, _edges[i].age);
		}
		return out;
	}

	/**
	 * Remove all edges older than min_age
	 */
	void purge_old_edges(float min_age) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].age > min_age)
				_edges.erase(_edges.begin()+i);
		}
	}

	/**
	 * Remove only the oldest edge, if older than min_age.
	 */
	void purge_oldest_edge(float min_age=0) {
		Lock ml(&mutex);
		unsigned idx = 0;
		XFLOAT idx_age = 0;
		for (unsigned i = 0; i < _edges.size(); i++) {
			float a = _edges[i].age;
			if (a >= idx_age)
			{
				idx = i;
				idx_age = a;
			}
		}
		if (idx_age > min_age)
			_edges.erase(_edges.begin()+idx);
	}

	/**
	 * Remove all edge-less nodes.
	 * Return them as a list.
	 */
	std::vector<unsigned> purge_orphans() {
		Lock ml(&mutex);
		std::vector<unsigned> orphans;
		for (unsigned ni = 0; ni < _nodes.size(); ni++) {
			if (_nodes[ni].active) {
				bool is_orphan = true;
				for (unsigned ei = 0; ei < _edges.size(); ei++)
					if (_edges[ei].n1 == ni || _edges[ei].n2 == ni) {
						is_orphan = false;
						break;
					}
				if (is_orphan)
					orphans.push_back(ni);
			}
		}

		for (unsigned i = 0; i < orphans.size(); i++)
			_remove_node(orphans[i]);

		return orphans;
	}

	/**
	 * Getters and setters.
	 */

	float get_edge_age(unsigned node1, unsigned node2) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node1 && _edges[i].n2 == node2 ||
			    _edges[i].n1 == node2 && _edges[i].n2 == node1)
				return _edges[i].age;
		}
		throw std::runtime_error("edge not found");
	}

	void set_edge_age(unsigned node1, unsigned node2, float age=0.) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++)
			if (_edges[i].n1 == node1 && _edges[i].n2 == node2 ||
			    _edges[i].n1 == node2 && _edges[i].n2 == node1) {
				_edges[i].age = age;
				return;
			}
		throw std::runtime_error("edge not found");
	}

	float get_node_error(unsigned node) {
		Lock ml(&mutex);
		return _nodes[node].error;
	}

	void set_node_error(unsigned node, float error=0.) {
		Lock ml(&mutex);
		_nodes[node].error = error;
	}

	void increment_node_error(unsigned node, float error) {
		Lock ml(&mutex);
		_nodes[node].error += error;
	}

	void increment_all_edge_ages(float age=1.) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++)
			_edges[i].age += age;
	}

	void increment_neighbour_edge_ages(unsigned node, float age=1.) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node || _edges[i].n2 == node)
				_edges[i].age += age;
		}
	}

	void increment_node_age(unsigned node, float amount=1) {
		Lock ml(&mutex);
		_nodes[node].age += amount;
	}

	float get_node_age(unsigned node) {
		Lock ml(&mutex);
		return _nodes[node].age;
	}

	void set_node_age(unsigned node, float age=0.) {
		Lock ml(&mutex);
		_nodes[node].age = age;
	}

	unsigned get_node_count() {
		Lock ml(&mutex);
		unsigned count = 0;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				count ++;
		return count;
	}

	unsigned get_edge_count() {
		Lock ml(&mutex);
		return _edges.size();
	}

	std::vector<unsigned> get_all_nodes() {
		Lock ml(&mutex);
		std::vector<unsigned> nodes;
		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				nodes.push_back(i);
		return nodes;
	}

	void reset_errors() {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _nodes.size(); i++)
			_nodes[i].error = 0.;
	}

	/**
	 * Return the oldest node.
	 */
	unsigned find_oldest_node() {
		Lock ml(&mutex);

		bool set = false;
		unsigned wpu;
		XFLOAT wpu_error = 0;
		for (unsigned i = 0; i < _nodes.size(); i++) {
			if (_nodes[i].active) {
				float e = _nodes[i].age;
				if (e >= wpu_error) {
					wpu = i;
					wpu_error = e;
					set = true;
				}
			}
		}

		if (!set)
			throw std::runtime_error("no node found");

		return wpu;
	}

	/**
	 * Return the worst performing unit.
	 */
	unsigned find_wpu() {
		Lock ml(&mutex);

		bool set = false;
		unsigned wpu;
		XFLOAT wpu_error = 0;
		for (unsigned i = 0; i < _nodes.size(); i++) {
			if (_nodes[i].active) {
				float e = _nodes[i].error;
				if (e >= wpu_error) {
					wpu = i;
					wpu_error = e;
					set = true;
				}
			}
		}

		if (!set)
			throw std::runtime_error("no wpu found");

		return wpu;
	}

	/**
	 * Return the worst performing neighbour of given node.
	 */
	unsigned find_wpu(unsigned node) {
		Lock ml(&mutex);
		std::vector<unsigned> neighbours = _get_neighbours(node);

		bool set = false;
		unsigned wpu;
		XFLOAT max_e = 0;
		for (unsigned i = 0; i < neighbours.size(); i++) {
			float e = _nodes[neighbours[i]].error;
			if (e >= max_e) {
				wpu = neighbours[i];
				max_e = e;
				set = true;
			}
		}

		if (!set)
			throw std::runtime_error("no neighbour wpu found");

		return wpu;
	}

	/**
	 * Return unique nodes connected to edge.
	 * Exclude node1 and node2.
	 */
	std::vector<unsigned> get_neighbourhood_of_edge(unsigned node1, unsigned node2)
	{
		std::vector<unsigned> n1 = _get_neighbours(node1);
		std::vector<unsigned> n2 = _get_neighbours(node2);
		std::vector<unsigned> n;

		for (unsigned i = 0; i < n1.size(); i++)
			if (n1[i] != node1 && n1[i] != node2 &&
			    std::find(n.begin(), n.end(), n1[i]) == n.end() )
				n.push_back(n1[i]);

		for (unsigned i = 0; i < n2.size(); i++)
			if (n2[i] != node1 && n2[i] != node2 &&
			    std::find(n.begin(), n.end(), n2[i]) == n.end() )
				n.push_back(n2[i]);

		return n;
	}

	void print_to_file(FileName &fn) {
		FILE * fp;
		fp = fopen (fn.c_str(), "w+");

		fprintf(fp, "_nodes [index error age]\n");

		for (unsigned i = 0; i < _nodes.size(); i++)
			if (_nodes[i].active)
				fprintf(fp, "%3d %5.1f %5.1f\n", i, _nodes[i].error, _nodes[i].age);

		fprintf(fp, "\n_edges [node1 node2 age]\n");

		for(unsigned i = 0; i < _edges.size(); i++)
			fprintf(fp, "%3d %3d %5.1f\n", _edges[i].n1, _edges[i].n2, _edges[i].age);

		fclose(fp);
	}
};

#endif
