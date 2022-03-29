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
 * Articlerank implementation. Uses the basic vertex-based API for
 * demonstration purposes.
 */

#include <string>
#include <fstream>
#include <cmath>

#define GRAPHCHI_DISABLE_COMPRESSION

#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;
 
#define RANDOMRESETPROB 0.15

typedef float VertexDataType;
typedef float EdgeDataType;

struct ArticleRankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
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
        float sum = 0;
        float edge_per_vertex = ((float)ginfo.nedges / (float)ginfo.nvertices);
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                edge->set_data(1.0 / v.num_outedges());
            }
            v.set_data(RANDOMRESETPROB); 
        } else {
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
            for(int i=0; i < v.num_inedges(); i++) {
                float val = v.inedge(i)->get_data();
                sum += val;                    
            }
            
            /* Compute my articlerank */
            float articlerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            /* Write my articlerank divided by the number of out-edges to
               each of my out-edges. */
            if (v.num_outedges() > 0) {
                float rank = articlerank / (v.num_outedges() + edge_per_vertex);
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
                    edge->set_data(rank);
                }
            }
            
            /* Set my new articlerank as the vertex value */
            v.set_data(articlerank); 
        }
    }  
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("articlerank");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 4);
    bool scheduler          = false;                    // Non-dynamic version of articlerank.
    int ntop                = get_option_int("top", 20);
    
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));

    /* Run */
    graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    
    ArticleRankProgram program;
    engine.run(program, niters);
    
    /* Output top ranked vertices */
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
    return 0;
}

