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

class SomGraph {

private:

	/**
	 * Class for graph nodes
	 */
	struct Node {
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

	std::unordered_map<unsigned, Node> _nodes;
	std::vector<Edge> _edges;

	pthread_mutex_t mutex;

	/*
	 * Non-thread-safe remove node.
	 */
	void _remove_node(unsigned node) {
		if (_nodes.find(node) == _nodes.end())
			throw std::runtime_error("node missing");

		for (unsigned i = 0; i < _edges.size(); i++)
			if (_edges[i].n1 == node || _edges[i].n2 == node)
				_edges.erase(_edges.begin()+i);

		_nodes.erase(node);
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
	 * Add an edge-less node to the graph.
	 */
	unsigned add_node() {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _nodes.size() + 1; i++)
			if (_nodes.find(i) == _nodes.end()) { // If index not found
				_nodes.emplace(i, Node{});
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
			throw std::runtime_error("cannot add edge to the node itself");

		if (_nodes.find(node1) == _nodes.end() || _nodes.find(node2) == _nodes.end())
			throw std::runtime_error("node missing");

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
		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n) {
			bool is_orphan = true;
			for (unsigned i = 0; i < _edges.size(); i++)
				if (_edges[i].n1 == n->first || _edges[i].n2 == n->first) {
					is_orphan = false;
					break;
				}
			if (is_orphan)
				orphans.push_back(n->first);
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

	void increment_neighbour_edge_ages(unsigned node) {
		Lock ml(&mutex);
		for (unsigned i = 0; i < _edges.size(); i++) {
			if (_edges[i].n1 == node || _edges[i].n2 == node)
				_edges[i].age ++;
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
		return _nodes.size();
	}

	unsigned get_edge_count() {
		Lock ml(&mutex);
		return _edges.size();
	}

	std::vector<unsigned> get_all_nodes() {
		Lock ml(&mutex);
		std::vector<unsigned> nodes;
		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n)
			nodes.push_back(n->first);
		return nodes;
	}

	void reset_errors() {
		Lock ml(&mutex);
		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n)
			n->second.error = 0.;
	}

	/**
	 * Return the oldest node.
	 */
	unsigned find_oldest_node() {
		Lock ml(&mutex);

		bool set = false;
		unsigned wpu;
		XFLOAT wpu_error = 0;
		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n) {
			float e = n->second.age;
			if (e >= wpu_error) {
				wpu = n->first;
				wpu_error = e;
				set = true;
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
		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n) {
			float e = n->second.error;
			if (e >= wpu_error) {
				wpu = n->first;
				wpu_error = e;
				set = true;
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
		std::ofstream out;
		out.open(fn);

		out << "_nodes [index error age]" << std::endl;

		for(std::unordered_map<unsigned, Node>::iterator n = _nodes.begin(); n != _nodes.end(); ++n)
			out << n->first << " " << n->second.error << " " << n->second.age << std::endl;

		out << std::endl;
		out << "_edges [node1 node2 age]" << std::endl;

		for(unsigned i = 0; i < _edges.size(); i++)
			out << _edges[i].n1 << " " << _edges[i].n2 << " " << _edges[i].age << std::endl;

		out.close();
	}
};

#endif
