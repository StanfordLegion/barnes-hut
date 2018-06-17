/* Copyright 2018 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "barnes_hut.h"

#include <stdlib.h>
#include "mappers/default_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

static LegionRuntime::Logger::Category log_mapper("barnes_hut");

class BarnesHutMapper : public DefaultMapper {
public:
  BarnesHutMapper(MapperRuntime *rt, Machine machine, Processor local, std::vector<Processor> *cpus, int sector_size);

  void slice_task(const MapperContext ctx,
                  const Task &task,
                  const SliceTaskInput &input,
                  SliceTaskOutput &output);

private:
  std::vector<Processor> &cpus;
  int sector_size;
};

BarnesHutMapper::BarnesHutMapper(MapperRuntime *rt, Machine machine, Processor local, std::vector<Processor> *cpus, int sector_size)
    : DefaultMapper(rt, machine, local, "barnes_hut_mapper"), cpus(*cpus), sector_size(sector_size) {
  std::set<Processor> all_procs;
  machine.get_all_processors(all_procs);
  // Recall that we create one mapper for every processor.  We
  // only want to print out this information one time, so only
  // do it if we are the mapper for the first processor in the
  // list of all processors in the machine.
  if (all_procs.begin()->id + 1 == local_proc.id) {
    // Print out how many processors there are and each
    // of their kinds.
    printf("There are %zd processors:\n", all_procs.size());
    for (std::set<Processor>::const_iterator it = all_procs.begin();
         it != all_procs.end(); it++) {
      // For every processor there is an associated kind
      Processor::Kind kind = it->kind();
      switch (kind) {
        // Latency-optimized cores (LOCs) are CPUs
        case Processor::LOC_PROC: {
          printf("  Processor ID " IDFMT " is CPU\n", it->id);
          break;
        }
          // Throughput-optimized cores (TOCs) are GPUs
        case Processor::TOC_PROC: {
          printf("  Processor ID " IDFMT " is GPU\n", it->id);
          break;
        }
          // Processor for doing I/O
        case Processor::IO_PROC: {
          printf("  Processor ID " IDFMT " is I/O Proc\n", it->id);
          break;
        }
          // Utility processors are helper processors for
          // running Legion runtime meta-level tasks and
          // should not be used for running application tasks
        case Processor::UTIL_PROC: {
          printf("  Processor ID " IDFMT " is utility\n", it->id);
          break;
        }
        default:
          assert(false);
      }
    }
    // We can also get the list of all the memories available
    // on the target architecture and print out their info.
    std::set<Memory> all_mems;
    machine.get_all_memories(all_mems);
    printf("There are %zd memories:\n", all_mems.size());
    for (std::set<Memory>::const_iterator it = all_mems.begin();
         it != all_mems.end(); it++) {
      Memory::Kind kind = it->kind();
      size_t memory_size_in_kb = it->capacity() >> 10;
      switch (kind) {
        // RDMA addressable memory when running with GASNet
        case Memory::GLOBAL_MEM: {
          printf("  GASNet Global Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // DRAM on a single node
        case Memory::SYSTEM_MEM: {
          printf("  System Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Pinned memory on a single node
        case Memory::REGDMA_MEM: {
          printf("  Pinned Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // A memory associated with a single socket
        case Memory::SOCKET_MEM: {
          printf("  Socket Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Zero-copy memory betweeen CPU DRAM and
          // all GPUs on a single node
        case Memory::Z_COPY_MEM: {
          printf("  Zero-Copy Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // GPU framebuffer memory for a single GPU
        case Memory::GPU_FB_MEM: {
          printf("  GPU Frame Buffer Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Disk memory on a single node
        case Memory::DISK_MEM: {
          printf("  Disk Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // HDF framebuffer memory for a single GPU
        case Memory::HDF_MEM: {
          printf("  HDF Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // File memory on a single node
        case Memory::FILE_MEM: {
          printf("  File Memory ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Block of memory sized for L3 cache
        case Memory::LEVEL3_CACHE: {
          printf("  Level 3 Cache ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Block of memory sized for L2 cache
        case Memory::LEVEL2_CACHE: {
          printf("  Level 2 Cache ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
          // Block of memory sized for L1 cache
        case Memory::LEVEL1_CACHE: {
          printf("  Level 1 Cache ID " IDFMT " has %zd KB\n",
                 it->id, memory_size_in_kb);
          break;
        }
        default:
          assert(false);
      }
    }

    // The Legion machine model represented by the machine object
    // can be thought of as a graph with processors and memories
    // as the two kinds of nodes.  There are two kinds of edges
    // in this graph: processor-memory edges and memory-memory
    // edges.  An edge between a processor and a memory indicates
    // that the processor can directly perform load and store
    // operations to that memory.  Memory-memory edges indicate
    // that data movement can be directly performed between the
    // two memories.  To illustrate how this works we examine
    // all the memories visible to our local processor in
    // this mapper.  We can get our set of visible memories
    // using the 'get_visible_memories' method on the machine.
    std::set<Memory> vis_mems;
    machine.get_visible_memories(local_proc, vis_mems);
    printf("There are %zd memories visible from processor " IDFMT "\n",
           vis_mems.size(), local_proc.id);
    for (std::set<Memory>::const_iterator it = vis_mems.begin();
         it != vis_mems.end(); it++) {
      // Edges between nodes are called affinities in the
      // machine model.  Affinities also come with approximate
      // indications of the latency and bandwidth between the
      // two nodes.  Right now these are unit-less measurements,
      // but our plan is to teach the Legion runtime to profile
      // these values on start-up to give them real values
      // and further increase the portability of Legion applications.
      std::vector<ProcessorMemoryAffinity> affinities;
      int results =
          machine.get_proc_mem_affinity(affinities, local_proc, *it);
      // We should only have found 1 results since we
      // explicitly specified both values.
      assert(results == 1);
      printf("  Memory " IDFMT " has bandwidth %d and latency %d\n",
             it->id, affinities[0].bandwidth, affinities[0].latency);
    }
  }
}

void BarnesHutMapper::slice_task(const MapperContext ctx,
                                 const Task &task,
                                 const SliceTaskInput &input,
                                 SliceTaskOutput &output) {
  const char *task_name = task.get_task_name();
  if (strcmp(task_name, "assign_sectors") != 0
      && strcmp(task_name, "update_body_force_root")
      && strcmp(task_name, "update_body_force") != 0
      && strcmp(task_name, "update_body_speed") != 0) {
    DefaultMapper::slice_task(ctx, task, input, output);
    return;
  }

//  printf("Task name %s\n", task_name);
//  printf("Volume %lu\n", input.domain.get_volume());

  double scaling_factor = (double) cpus.size() / sector_size;
//  printf("Scaling Factor %f\n", scaling_factor);

  output.slices.resize(input.domain.get_volume());
  unsigned idx = 0;
  for (Domain::DomainPointIterator itr(task.index_domain); itr; itr++, idx++) {
    TaskSlice &slice = output.slices[idx];
    slice.domain = Domain::from_point<1>(itr.p.get_point<1>());
    long sector = itr.p.get_point<1>().x[0];

//    printf("sector %lu\n", sector);
//    printf("cpu %f\n", sector * scaling_factor);

    slice.proc = cpus.at((int) sector * scaling_factor);
    slice.recurse = false;
    slice.stealable = false;
  }
}

static void create_mappers(Machine machine, HighLevelRuntime *runtime,
                           const std::set<Processor> &local_procs) {
  printf("Address space count %lu\n", machine.get_address_space_count());

  std::vector<Processor>* cpus = new std::vector<Processor>();;
  std::set<Processor> all_procs;
  machine.get_all_processors(all_procs);
  for (std::set<Processor>::const_iterator it = all_procs.begin();
       it != all_procs.end(); it++) {
    Processor::Kind kind = it->kind();
    switch (kind) {
      case Processor::LOC_PROC: {
        cpus->push_back(*it);
        break;
      }
      default:
        break;
    }
  }

  char* num_sectors_env = getenv("NUM_SECTORS");
  int num_sectors;
  if (num_sectors_env) {
    num_sectors = atoi(num_sectors_env);
  } else {
    num_sectors = 256;
  }

  for (std::set<Processor>::const_iterator it = local_procs.begin();
       it != local_procs.end(); it++) {
    BarnesHutMapper *mapper = new BarnesHutMapper(runtime->get_mapper_runtime(), machine, *it, cpus, num_sectors);
    runtime->replace_default_mapper(mapper, *it);
  }
}

void register_mappers() {
  HighLevelRuntime::add_registration_callback(create_mappers);
}



