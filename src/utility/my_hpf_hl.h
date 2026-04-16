#ifndef FLOW_MY_HPF_HL_H__
#define FLOW_MY_HPF_HL_H__

/* LICENSE
 *
 * The source code is subject to the following academic license.
 * Note this is not an open source license.
 *
 * Copyright © 2001. The Regents of the University of California (Regents).
 * All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * two paragraphs appear in all copies, modifications, and distributions.
 * Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck
 * Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial
 * licensing opportunities. Created by Bala Chandran and Dorit S. Hochbaum,
 * Department of Industrial Engineering and Operations Research,
 * University of California, Berkeley.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

#include <cstdint>
#include <cstdlib>

#include "memory_stack.h"
#include "macros.h"

namespace my_reimpls_hpf {
    enum TermType : uint32_t {
        SOURCE = 0,
        SINK = 1
    };

    inline uint64_t align64(uint64_t n) { return (n + 63) & ~uint64_t(63); }

    template<class T>
    class Array {
        T *data_ = nullptr;

        static T *assume_aligned(const T *p) {
            return static_cast<T *>(__builtin_assume_aligned(p, 64));
        }

    public:
        Array() = default;

        void init(const T *mem) { data_ = assume_aligned(mem); }
        T &operator[](size_t i) { return data_[i]; }
        const T &operator[](size_t i) const { return data_[i]; }
        T *data() { return data_; }
        const T *data() const { return data_; }
    };

    template<class Cap>
    struct HpfNode;

    template<class Cap>
    struct HpfArc {
        uint32_t from;
        uint32_t to;
        Cap flow;
        Cap capacity;
        bool direction;
    };

    template<class Cap>
    struct HpfNode {
        uint32_t numAdjacent;
        uint32_t label;
        Cap excess;

        HpfNode<Cap> *parent;
        HpfNode<Cap> *childList;
        HpfNode<Cap> *nextScan;

        uint32_t numOutOfTree;
        HpfArc<Cap> **outOfTree;
        uint32_t nextArc;
        HpfArc<Cap> *arcToParent;
        HpfNode<Cap> *next;

        void addOutOfTree(HpfArc<Cap> *out);
    };

    template<class Cap>
    struct HpfRoot {
        HpfNode<Cap> *start;
    };

    template<class Cap>
    class Hpf {
    public:
        using Node = HpfNode<Cap>;
        using Arc = HpfArc<Cap>;
        using Root = HpfRoot<Cap>;

        static constexpr TermType source_side = SOURCE;
        static constexpr TermType sink_side = SINK;

        Array<Arc> arcList;

        Hpf(size_t expectedNodes, size_t expectedArcs, HeiProMap::MemoryStack &mem_stack);

        uint32_t add_node(uint32_t num = 1);

        void add_edge(uint32_t from, uint32_t to, Cap capacity);

        void mincut();

        TermType what_label(uint32_t node) const;

        Cap compute_maxflow() const noexcept;

        uint32_t get_num_arcs() const noexcept { return numArcs; }

        inline void set_source(uint32_t s) { source = s; }
        inline void set_sink(uint32_t t) { sink = t; }

        void reset(size_t expectedNodes, size_t expectedArcs, HeiProMap::MemoryStack &mem_stack) {
            numNodes = 0;
            numArcs = 0;
            source = 0;
            sink = 0;
            highestStrongLabel = 1;

            const uint64_t node_mem = align64(sizeof(Node) * expectedNodes);
            const uint64_t root_mem = align64(sizeof(Root) * expectedNodes);
            const uint64_t label_mem = align64(sizeof(uint32_t) * expectedNodes);
            const uint64_t arc_mem = align64(sizeof(Arc) * expectedArcs);
            const uint64_t out_mem = align64(sizeof(Arc *) * expectedArcs * 2ULL);
            mem_stack.ensure(node_mem + root_mem + label_mem + arc_mem + out_mem);
            mem_stack.clear();

            adjacencyList.init(static_cast<Node *>(mem_stack.get_memory(node_mem)));
            strongRoots.init(static_cast<Root *>(mem_stack.get_memory(root_mem)));
            labelCount.init(static_cast<uint32_t *>(mem_stack.get_memory(label_mem)));
            arcList.init(static_cast<Arc *>(mem_stack.get_memory(arc_mem)));
            outOfTreePtrs.init(static_cast<Arc **>(mem_stack.get_memory(out_mem)));
        }

    private:
        uint32_t numNodes;
        uint32_t numArcs;
        uint32_t source;
        uint32_t sink;

        uint32_t highestStrongLabel;

        Array<Node> adjacencyList;
        Array<Root> strongRoots;
        Array<uint32_t> labelCount;
        Array<Arc *> outOfTreePtrs;

        void init_mincut();

        void addToStrongBucket(Node *newRoot, Root *rootBucket);

        Node *getHighestStrongRoot();

        Node *getNextStrongRoot();

        void processRoot(Node *strongRoot);

        Arc *findWeakNode(Node *strongNode, Node **weakNode);

        void merge(Node *parent, Node *child, Arc *newArc);

        void addRelationship(Node *newParent, Node *child);

        void breakRelationship(Node *oldParent, Node *child);

        void pushExcess(Node *strongRoot);

        void pushUpward(Arc *currentArc, Node *child, Node *parent, Cap resCap);

        void pushDownward(Arc *currentArc, Node *child, Node *parent, Cap flow);

        void checkChildren(Node *curNode);

        void liftAll(Node *rootNode);
    };

    template<class Cap>
    inline Hpf<Cap>::Hpf(size_t expectedNodes, size_t expectedArcs, HeiProMap::MemoryStack &mem_stack) : arcList(),
                                                                                              numNodes(0),
                                                                                              numArcs(0),
                                                                                              source(0),
                                                                                              sink(0),
                                                                                              highestStrongLabel(1),
                                                                                              adjacencyList(),
                                                                                              strongRoots(),
                                                                                              labelCount(),
                                                                                              outOfTreePtrs() {
        reset(expectedNodes, expectedArcs, mem_stack);
    }

    template<class Cap>
    inline uint32_t Hpf<Cap>::add_node(uint32_t num) {
        numNodes += num;
        for (uint32_t i = 0; i < numNodes; ++i) {
            adjacencyList[i] = Node{};
            strongRoots[i] = Root{};
            labelCount[i] = 0;
        }
        return numNodes;
    }

    template<class Cap>
    inline void Hpf<Cap>::add_edge(uint32_t from, uint32_t to, Cap capacity) {
        Arc &arc = arcList[numArcs];
        arc.from = from;
        arc.to = to;
        arc.flow = 0;
        arc.capacity = capacity;
        arc.direction = true;
        adjacencyList[from].numAdjacent += 1;
        adjacencyList[to].numAdjacent += 1;
        numArcs++;
    }

    template<class Cap>
    inline void Hpf<Cap>::mincut() {
        init_mincut();

        // pseudoflowPhase1
        Node *strongRoot;

        while ((strongRoot = getNextStrongRoot())) {
            processRoot(strongRoot);
        }
    }

    template<class Cap>
    inline TermType Hpf<Cap>::what_label(uint32_t node) const {
        return adjacencyList[node].label >= numNodes ? SOURCE : SINK;
    }

    template<class Cap>
    inline Cap Hpf<Cap>::compute_maxflow() const noexcept {
        Cap cut = 0;

        // Compute value of minimum cut which is equal to the max. flow
        for (uint32_t i = 0; i < numArcs; ++i) {
            const Arc &a = arcList[i];
            if (adjacencyList[a.from].label >= numNodes && adjacencyList[a.to].label < numNodes) {
                cut += a.capacity;
            }
        }

        return cut;
    }

    template<class Cap>
    inline void Hpf<Cap>::init_mincut() {
        Arc **crntOutOfTree = outOfTreePtrs.data();
        for (uint32_t i = 0; i < numNodes; ++i) {
            Node &n = adjacencyList[i];
            // createOutOfTree
            if (n.numAdjacent) {
                n.outOfTree = crntOutOfTree;
                crntOutOfTree += n.numAdjacent;
            }
        }

        for (uint32_t i = 0; i < numArcs; ++i) {
            Arc &a = arcList[i];
            uint32_t from = a.from;
            uint32_t to = a.to;
            if (!(source == to || sink == from || from == to)) {
                if (source == from && to == sink) {
                    a.flow = a.capacity;
                } else if (from == source) {
                    adjacencyList[from].addOutOfTree(&a);
                } else if (to == sink) {
                    adjacencyList[to].addOutOfTree(&a);
                } else {
                    adjacencyList[from].addOutOfTree(&a);
                }
            }
        }

        // simpleInitialization
        for (uint32_t i = 0; i < adjacencyList[source].numOutOfTree; ++i) {
            Arc *tempArc = adjacencyList[source].outOfTree[i];
            tempArc->flow = tempArc->capacity;
            adjacencyList[tempArc->to].excess += tempArc->capacity;
        }

        for (uint32_t i = 0; i < adjacencyList[sink].numOutOfTree; ++i) {
            Arc *tempArc = adjacencyList[sink].outOfTree[i];
            tempArc->flow = tempArc->capacity;
            adjacencyList[tempArc->from].excess -= tempArc->capacity;
        }

        adjacencyList[source].excess = 0;
        adjacencyList[sink].excess = 0;

        for (uint32_t i = 0; i < numNodes; ++i) {
            if (adjacencyList[i].excess > 0) {
                adjacencyList[i].label = 1;
                ++labelCount[1];

                addToStrongBucket(&adjacencyList[i], &strongRoots[1]);
            }
        }

        adjacencyList[source].label = numNodes;
        adjacencyList[sink].label = 0;
        labelCount[0] = (numNodes - 2) - labelCount[1];
    }

    template<class Cap>
    inline void Hpf<Cap>::addToStrongBucket(Node *newRoot, Root *rootBucket) {
        newRoot->next = rootBucket->start;
        rootBucket->start = newRoot;
    }

    template<class Cap>
    inline typename Hpf<Cap>::Node *
    Hpf<Cap>::getHighestStrongRoot() {
        Node *strongRoot;

        for (uint32_t i = highestStrongLabel; i > 0; --i) {
            if (strongRoots[i].start) {
                highestStrongLabel = i;
                if (labelCount[i - 1]) {
                    strongRoot = strongRoots[i].start;
                    strongRoots[i].start = strongRoot->next;
                    strongRoot->next = nullptr;
                    return strongRoot;
                }

                while (strongRoots[i].start) {
                    strongRoot = strongRoots[i].start;
                    strongRoots[i].start = strongRoot->next;
                    liftAll(strongRoot);
                }
            }
        }

        if (!strongRoots[0].start) {
            return nullptr;
        }

        while (strongRoots[0].start) {
            strongRoot = strongRoots[0].start;
            strongRoots[0].start = strongRoot->next;
            strongRoot->label = 1;
            --labelCount[0];
            ++labelCount[1];

            addToStrongBucket(strongRoot, &strongRoots[strongRoot->label]);
        }

        highestStrongLabel = 1;

        strongRoot = strongRoots[1].start;
        strongRoots[1].start = strongRoot->next;
        strongRoot->next = nullptr;

        return strongRoot;
    }

    template<class Cap>
    inline typename Hpf<Cap>::Node *
    Hpf<Cap>::getNextStrongRoot() {
        return getHighestStrongRoot();
    }

    template<class Cap>
    inline void Hpf<Cap>::processRoot(Node *strongRoot) {
        Node *strongNode = strongRoot, *weakNode;
        Arc *out;

        strongRoot->nextScan = strongRoot->childList;

        if ((out = findWeakNode(strongRoot, &weakNode))) {
            merge(weakNode, strongNode, out);
            pushExcess(strongRoot);
            return;
        }

        checkChildren(strongRoot);

        while (strongNode) {
            while (strongNode->nextScan) {
                Node *temp = strongNode->nextScan;
                strongNode->nextScan = strongNode->nextScan->next;
                strongNode = temp;
                strongNode->nextScan = strongNode->childList;

                if ((out = findWeakNode(strongNode, &weakNode))) {
                    merge(weakNode, strongNode, out);
                    pushExcess(strongRoot);
                    return;
                }

                checkChildren(strongNode);
            }

            if ((strongNode = strongNode->parent)) {
                checkChildren(strongNode);
            }
        }

        addToStrongBucket(strongRoot, &strongRoots[strongRoot->label]);

        ++highestStrongLabel;
    }

    template<class Cap>
    inline typename Hpf<Cap>::Arc *Hpf<Cap>::findWeakNode(Node *strongNode, Node **weakNode) {
        Arc *out;

        const uint32_t extremumStrongLabel = highestStrongLabel;
        uint32_t size = strongNode->numOutOfTree;

        for (uint32_t i = strongNode->nextArc; i < size; ++i) {
            if (adjacencyList[strongNode->outOfTree[i]->to].label == (extremumStrongLabel - 1)) {
                strongNode->nextArc = i;
                out = strongNode->outOfTree[i];
                (*weakNode) = &adjacencyList[out->to];
                --strongNode->numOutOfTree;
                strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
                return out;
            }

            if (adjacencyList[strongNode->outOfTree[i]->from].label == (extremumStrongLabel - 1)) {
                strongNode->nextArc = i;
                out = strongNode->outOfTree[i];
                (*weakNode) = &adjacencyList[out->from];
                --strongNode->numOutOfTree;
                strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
                return out;
            }
        }

        strongNode->nextArc = strongNode->numOutOfTree;

        return nullptr;
    }

    template<class Cap>
    inline void Hpf<Cap>::merge(Node *parent, Node *child, Arc *newArc) {
        Node *current = child, *newParent = parent;

        while (current->parent) {
            Arc *oldArc = current->arcToParent;
            current->arcToParent = newArc;
            Node *oldParent = current->parent;
            breakRelationship(oldParent, current);
            addRelationship(newParent, current);
            newParent = current;
            current = oldParent;
            newArc = oldArc;
            newArc->direction = 1 - newArc->direction;
        }

        current->arcToParent = newArc;
        addRelationship(newParent, current);
    }

    template<class Cap>
    inline void Hpf<Cap>::addRelationship(Node *newParent, Node *child) {
        child->parent = newParent;
        child->next = newParent->childList;
        newParent->childList = child;
    }

    template<class Cap>
    inline void Hpf<Cap>::breakRelationship(Node *oldParent, Node *child) {
        Node *current;

        child->parent = nullptr;
        if (oldParent->childList == child) {
            oldParent->childList = child->next;
            child->next = nullptr;
            return;
        }

        for (current = oldParent->childList; current->next != child; current = current->next) {
            // Do nothing
        }

        current->next = child->next;
        child->next = nullptr;
    }

    template<class Cap>
    inline void Hpf<Cap>::pushExcess(Node *strongRoot) {
        Node *current, *parent;
        int prevEx = 1;

        for (current = strongRoot; current->excess && current->parent; current = parent) {
            parent = current->parent;
            prevEx = parent->excess;

            Arc *arcToParent = current->arcToParent;

            if (arcToParent->direction) {
                pushUpward(arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow));
            } else {
                pushDownward(arcToParent, current, parent, arcToParent->flow);
            }
        }

        if ((current->excess > 0) && (prevEx <= 0)) {
            addToStrongBucket(current, &strongRoots[current->label]);
        }
    }

    template<class Cap>
    inline void Hpf<Cap>::pushUpward(Arc *currentArc, Node *child, Node *parent, Cap resCap) {
        if (resCap >= child->excess) {
            parent->excess += child->excess;
            currentArc->flow += child->excess;
            child->excess = 0;
            return;
        }

        currentArc->direction = 0;
        parent->excess += resCap;
        child->excess -= resCap;
        currentArc->flow = currentArc->capacity;
        parent->outOfTree[parent->numOutOfTree] = currentArc;
        ++parent->numOutOfTree;
        breakRelationship(parent, child);

        addToStrongBucket(child, &strongRoots[child->label]);
    }

    template<class Cap>
    inline void Hpf<Cap>::pushDownward(Arc *currentArc, Node *child, Node *parent, Cap flow) {
        if (flow >= child->excess) {
            parent->excess += child->excess;
            currentArc->flow -= child->excess;
            child->excess = 0;
            return;
        }

        currentArc->direction = 1;
        child->excess -= flow;
        parent->excess += flow;
        currentArc->flow = 0;
        parent->outOfTree[parent->numOutOfTree] = currentArc;
        ++parent->numOutOfTree;
        breakRelationship(parent, child);

        addToStrongBucket(child, &strongRoots[child->label]);
    }

    template<class Cap>
    inline void Hpf<Cap>::checkChildren(Node *curNode) {
        for (; curNode->nextScan; curNode->nextScan = curNode->nextScan->next) {
            if (curNode->nextScan->label == curNode->label) {
                return;
            }
        }

        --labelCount[curNode->label];
        ++curNode->label;
        ++labelCount[curNode->label];

        curNode->nextArc = 0;
    }

    template<class Cap>
    inline void Hpf<Cap>::liftAll(Node *rootNode) {
        Node *current = rootNode;

        current->nextScan = current->childList;

        --labelCount[current->label];
        current->label = numNodes;

        for (; (current); current = current->parent) {
            while (current->nextScan) {
                Node *temp = current->nextScan;
                current->nextScan = current->nextScan->next;
                current = temp;
                current->nextScan = current->childList;

                --labelCount[current->label];
                current->label = numNodes;
            }
        }
    }

    template<class Cap>
    inline void HpfNode<Cap>::addOutOfTree(HpfArc<Cap> *out) {
        outOfTree[numOutOfTree] = out;
        numOutOfTree++;
    }
} // namespace my_reimpls_hpf

#endif // FLOW_MY_HPF_HL_H__
