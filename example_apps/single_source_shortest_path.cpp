
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
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */



#include <string>
#include <float.h>
#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;
/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
typedef float VertexDataType;
typedef float EdgeDataType;

/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */
struct SSSPProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        
        if (gcontext.iteration == 0) {
            /* Vertex init */
            if(vertex.id() == 0){
                vertex.set_data(0);
            }else{
                vertex.set_data(FLT_MAX);
            }
            /* Outedge init */
            for(int i=0; i < vertex.num_outedges(); i++) {
                graphchi_edge<float> * edge = vertex.outedge(i);
                edge->set_data(FLT_MAX);
            }
        } else {
            /* Loop over in-edges */
            float curr_min_dist = vertex.get_data();
            for(int i=0; i < vertex.num_inedges(); i++) {
                float min_dist = vertex.inedge(i)->get_data();
                if(min_dist != FLT_MAX && min_dist < curr_min_dist){
                    curr_min_dist = min_dist;
                }
            }
            
            /* Loop over out-edges */
            for(int i=0; i < vertex.num_outedges(); i++) {
                vertex.outedge(i)->set_data(curr_min_dist+1);
            }
            vertex.set_data(curr_min_dist);
        }
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
    }
    
    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext) {
    }
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("sssp");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 4); // Number of iterations
    bool scheduler       = false; // Whether to use selective scheduling
    int ntop             = get_option_int("top", 20);

    /* Detect the number of shards or preprocess an input to create them */
    int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
    
    /* Run */
    SSSPProgram program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.run(program, niters);
    
    /* Output top ranked vertices */
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
