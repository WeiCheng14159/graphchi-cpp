/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>

#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;

typedef uint64_t VertexDataType;
typedef uint64_t EdgeDataType;


struct TriangleCountingProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {

    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &ginfo) {
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
        if (iteration == 4) {
            std::cout << "Converged!" << std::endl;
            ginfo.set_last_iteration(iteration);
        }
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
            
        if (ginfo.iteration == 0) { // Clear all edge data
            for(int i=0; i < v.num_edges(); i++) {
                graphchi_edge<uint64_t> * edge = v.edge(i);
                edge->set_data(ginfo.nvertices + 1);
            }
            v.set_data(0);
        } else if(ginfo.iteration == 1 ) { // Hop 1 first step
            for(int i=0; i < v.num_edges(); i++) {
                graphchi_edge<uint64_t> * edge = v.edge(i);
                // Sends its ID to all neighbors with higher ID than its own.
                if(edge->vertex_id() > v.id()){
                    uint64_t new_edata = 1 << (v.id() % 64);
                    edge->set_data(new_edata);
                    std::cout << "[Iter 1]" << v.id() << "->" << edge->vertex_id() << ":" << new_edata << std::endl;
                }
            }
        } else if(ginfo.iteration == 2) { // Hop 1 second step
            uint64_t in_edata = 0;
            // Receive from lower vertices
            for(int i=0; i < v.num_edges(); i++) {
                graphchi_edge<uint64_t> * edge = v.edge(i);
                if(edge->vertex_id() < v.id()) {
                     in_edata |= edge->get_data();
                }
            }
            // Add my vertex ID
            uint64_t new_v = (1 << (v.id() % 64)) | in_edata;
            v.set_data(new_v);
            std::cout << "[Iter 2]" << v.id() << " get new v = " << new_v << std::endl;
        } else if(ginfo.iteration == 3) { // Hop 2 first step
            for(int i=0; i < v.num_edges(); i++) {
                graphchi_edge<uint64_t> * edge = v.edge(i);
                // Sends its vertex weights to all neighbors with higher ID than its own.
                if(edge->vertex_id() > v.id()){
                    edge->set_data(v.get_data());
                    std::cout << "[Iter 3]" << v.id() << "->" << edge->vertex_id() << ":" << v.get_data() << std::endl;
                }
            }
        } else if(ginfo.iteration == 4) { // Hop 2 second step
            uint64_t adj_hash = 0;
            for(int j=0; j < v.num_edges(); j++) {
                graphchi_edge<uint64_t> * edge = v.edge(j);
                if(edge->vertex_id() < v.id()){
                    adj_hash |= (1 << (edge->vertex_id() % 64));
                }
            }
            std::cout << "[Iter 4]" << v.id() << ":adj_hash=" << adj_hash << std::endl;

            uint64_t in_edata = 0;
            for(int j=0; j < v.num_edges(); j++) {
                graphchi_edge<uint64_t> * edge = v.edge(j);
                if(edge->vertex_id() < v.id()){
                    in_edata |= edge->get_data();
                }
            }
            std::cout << "[Iter 4]"<< v.id() << ":in_edata=" << in_edata << std::endl;
            
            uint64_t and_op = adj_hash & in_edata;
            uint32_t pop_count =  __builtin_popcountl(and_op);
            VertexDataType tri_cnt = pop_count * (pop_count-1) / 2;
            v.set_data(tri_cnt);
            std::cout << v.id() << ":tri_cnt=" << tri_cnt << std::endl;
        }
    }
    
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("simple-triangle-counting");
    
    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 4);
    bool scheduler          = false;                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 20);
    
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));

    /* Run */
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    
    TriangleCountingProgram program;
    engine.run(program, niters);
    
    /* Output top ranked vertices */
    std::vector< vertex_value<VertexDataType> > top = get_top_vertices<VertexDataType>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
    return 0;
}

