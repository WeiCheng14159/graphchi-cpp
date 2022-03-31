
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
 * BFS algorithm implements in graphchi
 */



#include <string>
#include <float.h>
#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

#define START_VERTEX 0

using namespace graphchi;

typedef float VertexDataType;
typedef float EdgeDataType;

/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */
struct BFSProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    bool converged;
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        /* Remove itself from scheduler */
        gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());

        if (gcontext.iteration == 0) {
            /* Every vertex should run in the first iteration */
            gcontext.scheduler->add_task(vertex.id());
            /* Vertex init */
            if(vertex.id() == START_VERTEX){
                vertex.set_data(0); // Explored
                for(int i=0; i < vertex.num_outedges(); i++) { // Outedge init
                    graphchi_edge<float> * edge = vertex.outedge(i);
                    edge->set_data(1);
                }
            }else{
                vertex.set_data(FLT_MAX); // Not explored
                for(int i=0; i < vertex.num_outedges(); i++) { // Outedge init
                    graphchi_edge<float> * edge = vertex.outedge(i);
                    edge->set_data(FLT_MAX);
                }
            }
        } else {
            float curr_min_dist = vertex.get_data();
            /* Loop over in-edges */
            for(int i=0; i < vertex.num_inedges(); i++) {
                float min_dist = vertex.inedge(i)->get_data();
                if(min_dist < curr_min_dist){
                    curr_min_dist = min_dist;
                    converged = false;
                }
            }

            if(!converged){
                /* Loop over out-edges */
                for(int i=0; i < vertex.num_outedges(); i++) {
                    graphchi_edge<float> * edge = vertex.outedge(i);
                    edge->set_data(curr_min_dist+1);
                    gcontext.scheduler->add_task(edge->vertex_id());
                }
                vertex.set_data(curr_min_dist);
            }   
        }
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
        converged = iteration > 0;
    }
    
    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext) {
        if (converged) {
            std::cout << "Converged!" << iteration<<std::endl;
            gcontext.set_last_iteration(iteration);
        }
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
    metrics m("bfs");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 4); // Number of iterations
    bool scheduler       = true; // Whether to use selective scheduling
    int ntop             = get_option_int("top", 20);

    /* Detect the number of shards or preprocess an input to create them */
    int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
    
    /* Run */
    BFSProgram program;
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
