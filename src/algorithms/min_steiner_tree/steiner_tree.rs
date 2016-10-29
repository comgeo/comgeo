use traits::*;
use std::fmt;
use std::slice::{Iter, IterMut};
use std::iter::*;

#[derive(Debug)]
pub struct SteinerTree<P> {
    nodes: Vec<Node<P>>,
    terminal_count: usize
}

#[derive(Debug)]
pub struct Node<P> {
    p: P,
    ns: Vec<*mut Node<P>>,
    is_terminal: bool,
    id: usize
}

impl<P> Node<P> {
    fn new(p: P, is_terminal: bool) -> Self {
        Node {
            p: p,
            ns: Vec::new(),
            is_terminal: is_terminal,
            id: 0
        }
    }

    fn add_neighbour(&mut self, n: *mut Self) {
        self.ns.push(n);
    }

    fn remove_neighbour(&mut self, id: usize) {
        for i in 0..self.ns.len() {
            let neighbour = unsafe { &*self.ns[i] };
            if neighbour.id == id {
                self.ns.swap_remove(i);
                return;
            }
        }
    }


}

impl<P: Clone> Node<P> {
    fn clone(&self) -> Self {
        let mut ns = Vec::with_capacity(self.ns.len());
        for &n in self.ns.iter() {
            ns.push(n);
        }

        Node {
            p: self.p.clone(),
            ns: ns,
            is_terminal: self.is_terminal,
            id: self.id
        }
    }
}


pub struct Edges<'a, P: 'a> {
    nodes: Iter<'a, Node<P>>,
    edges: Option<IncidentEdges<'a, P>>
}
impl<'a, P> Iterator for Edges<'a, P> {
    type Item = Edge<'a, P>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(ref mut es) = self.edges {
                if let se @ Some(_) = es.filter(|e| e.nodes().0.id() < e.nodes().1.id()).next() {
                    return se
                }
            }
            match self.nodes.next().map(|n| n.edges()) {
                None => return None,
                a @ Some(_) => self.edges = a,
            }
        }
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.nodes.len() , Some(self.nodes.len()))
    }
}
impl<'a, P> ExactSizeIterator for Edges<'a, P> {}

pub struct Terminals<'a, P: 'a> {
    iter: Iter<'a, Node<P>>,
    terms: usize
}
impl<'a, P> Iterator for Terminals<'a, P> {
    type Item = &'a Node<P>;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(n) = self.iter.next() {
            if n.is_terminal() {
                return Some(n);
            } else {
                continue;
            }
        }

        None
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.terms, Some(self.terms))
    }
}
impl<'a, P> ExactSizeIterator for Terminals<'a, P> {}

pub struct SteinerPoints<'a, P: 'a> {
    iter: IterMut<'a, Node<P>>,
    steiner_points: usize
}
impl<'a, P> Iterator for SteinerPoints<'a, P> {
    type Item = &'a mut Node<P>;
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(n) = self.iter.next() {
            if !n.is_terminal() {
                return Some(n);
            } else {
                continue;
            }
        }

        None
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.steiner_points, Some(self.steiner_points))
    }
}
impl<'a, P> ExactSizeIterator for SteinerPoints<'a, P> {}


pub struct NodePairs<'a, P: 'a> {
    nodes: &'a [Node<P>],
    i1: usize,
    i2: usize
}
impl<'a, P> NodePairs<'a, P> {
    fn new(n: &'a [Node<P>]) -> NodePairs<'a, P> {
        NodePairs {
            nodes: n,
            i1: 0,
            i2: 1
        }
    }
}
impl<'a, P> Iterator for NodePairs<'a, P> {
    type Item = (&'a Node<P>, &'a Node<P>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i1 == self.nodes.len() - 1 {
            return None;
        }
        let ret = unsafe {
            (self.nodes.get_unchecked(self.i1),
             self.nodes.get_unchecked(self.i2))
        };
        self.i2 += 1;
        if self.i2 == self.nodes.len() {
            self.i1 += 1;
            self.i2 = self.i1 + 1;
        }
        Some(ret)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.nodes.len()*(self.nodes.len()-1)/2,
         Some(self.nodes.len()*(self.nodes.len()-1)/2))
    }
}
impl<'a, P> ExactSizeIterator for NodePairs<'a, P> {}

impl<P: Point> SteinerTree<P> {
    pub fn new(t: &[P], s: &[P], edges: &[(usize, usize)]) -> Self {
        let mut st = SteinerTree {
            nodes: t.iter().map(|x| Node::new(x.clone(), true)).chain(
                       s.iter().map(|x| Node::new(x.clone(), false))
                   ).collect(),
            terminal_count: t.len()
        };

        // TODO: Hack!! <- neighbour pointers are invalidated on reallocation
        st.nodes.reserve(100);

        for i in 0..st.nodes().len() {
            st.nodes[i].id = i;
        }

        for &(a, b) in edges.iter() {
            debug_assert!(a < t.len() + s.len());
            debug_assert!(b < t.len() + s.len());
            debug_assert!(a != b);

            let aptr = &mut st.nodes[a] as *mut Node<P>;
            let bptr = &mut st.nodes[b] as *mut Node<P>;
            st.nodes[a].add_neighbour(bptr);
            st.nodes[b].add_neighbour(aptr);
        }

        st
    }

    pub fn edges(&self) -> Edges<P> {
        Edges {
            nodes: self.nodes.iter(),
            edges: None
        }
    }

    pub fn node_pairs(&self) -> NodePairs<P> {
        NodePairs::new(&self.nodes[..])
    }

    pub fn nodes(&self) -> Iter<Node<P>> {
        self.nodes.iter()
    }

    pub fn terminals(&self) -> Terminals<P> {
        Terminals {
            iter: self.nodes.iter(),
            terms: self.terminal_count
        }
    }

    pub fn nodes_mut(&mut self) -> IterMut<Node<P>> {
        self.nodes.iter_mut()
    }

    pub fn steiner_points(&mut self) -> SteinerPoints<P> {
        SteinerPoints {
            steiner_points: self.nodes.len() - self.terminal_count,
            iter: self.nodes.iter_mut()
        }
    }

    pub fn i(&self, i: usize) -> &Node<P> {
        &self.nodes[i]
    }

    pub fn steiner_i(&mut self, i: usize) -> &mut Node<P> {
        //println!("steiner_i: i={}", i);
        &mut self.nodes[i]
    }

    pub fn non_degenerate(&mut self) {
        //TODO
    }

    pub fn remove_edge(&mut self, (a, b): (usize, usize)) {
        //println!("remove_edge: ({}, {})", a, b);
        //stdout().flush();
        debug_assert!(a < self.nodes.len());
        debug_assert!(b < self.nodes.len());

        self.nodes[a].remove_neighbour(b);
        self.nodes[b].remove_neighbour(a);
    }

    pub fn add_edge(&mut self, (a, b): (usize, usize)) {
        //println!("add_edge: ({}, {})", a, b);
        //stdout().flush();
        debug_assert!(a < self.nodes.len());
        debug_assert!(b < self.nodes.len());

        let aptr = &mut self.nodes[a] as *mut Node<P>;
        let bptr = &mut self.nodes[b] as *mut Node<P>;
        self.nodes[a].add_neighbour(bptr);
        self.nodes[b].add_neighbour(aptr);
    }

    pub fn push_node(&mut self, p: P, ns: &[usize], is_terminal: bool) -> &Node<P> {
        //println!("add_node: [ {:?} ]", ns);
        //stdout().flush();
        let mut node = Node::new(p, is_terminal);
        node.id = self.nodes.len();

        for &n in ns {
            debug_assert!(n < self.nodes.len());
            node.add_neighbour(&mut self.nodes[n]);
        }
        self.nodes.push(node);
        for &n in ns {
            let ptr = self.nodes.last_mut().unwrap() as *mut Node<P>;
            self.nodes[n].add_neighbour(ptr);
        }

        if is_terminal {
            self.terminal_count += 1;
        } else {
            self.nodes.last_mut().unwrap().init();
        }

        self.nodes.last().unwrap()
    }

    pub fn pop_node(&mut self) -> Option<P> {
        //println!("remove_node: {}", id);
        //stdout().flush();

        if let Some(node) = self.nodes.last() {
            if node.is_terminal() {
                self.terminal_count -= 1;
            }

            // Remove pointers to the node
            for &ptr in node.ns.iter() {
                let n = unsafe { &mut *ptr };
                n.remove_neighbour(node.id());
            }
        }

        self.nodes.pop().map(|node| node.p)
    }

    pub fn last_node(&self) -> Option<&Node<P>> {
        self.nodes.last()
    }
}

impl<P: Point> SteinerTree<P> {
    pub fn len<M: MinkowskiSpace<P>>(&self, geo: &M) -> P::R {
        self.edges().fold(P::R::zero(), |acc, e| acc + e.len(geo))
    }
}

impl<P: Point> Clone for SteinerTree<P> {
    fn clone(&self) -> Self {
        let mut nodes = Vec::with_capacity(self.nodes.len());
        for node in self.nodes() {
            nodes.push(node.clone());
        }

        for i in 0..nodes.len() {
            for n in 0..nodes[i].ns.len() {
                let id = unsafe { (*nodes[i].ns[n]).id() };
                let ptr = &mut nodes[id] as *mut Node<P>;
                nodes[i].ns[n] = ptr;
            }
        }

        SteinerTree {
            nodes: nodes,
            terminal_count: self.terminal_count
        }
    }
}

impl<E: Point + fmt::Display> fmt::Display for SteinerTree<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "Steiner tree:"));
        try!(writeln!(f, "\tTerminals:"));
        for n in self.terminals() {
            try!(writeln!(f, "\t  #{}: {}", n.id(), n));
        }
        try!(writeln!(f, "\tSteiner points:"));
        for n in self.nodes().filter(|n| !n.is_terminal()) {
            try!(writeln!(f, "\t  #{}: {}", n.id(), n));
        }
        Ok(())
    }
}





#[derive(Debug)]
pub struct Edge<'a, P: 'a> {
    n1: &'a Node<P>,
    n2: &'a Node<P>
}

impl<'a, P: 'a> Edge<'a, P> {
    pub fn nodes(&self) -> (&Node<P>, &Node<P>) {
        (self.n1, self.n2)
    }
}

impl<'a, P: 'a + Point> Edge<'a, P> {
    pub fn len<M: MinkowskiSpace<P>>(&self, geo: &M) -> P::R {
        geo.dist(self.n1.p(), self.n2.p())
    }
}

pub struct NeighboursMut<'a, P: 'a> {
    iter: Iter<'a, *mut Node<P>>,
}
impl<'a, P> Iterator for NeighboursMut<'a, P> {
    type Item = &'a mut Node<P>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|&n| unsafe { &mut *n })
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}
impl<'a, P> ExactSizeIterator for NeighboursMut<'a, P> {}

pub struct Neighbours<'a, P: 'a> {
    iter: Iter<'a, *mut Node<P>>,
}
impl<'a, P> Iterator for Neighbours<'a, P> {
    type Item = &'a Node<P>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|&n| unsafe { &*n })
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}
impl<'a, P> ExactSizeIterator for Neighbours<'a, P> {}

pub struct NeighboursData<'a, P: 'a> {
    iter: Iter<'a, *mut Node<P>>,
}
impl<'a, P> Iterator for NeighboursData<'a, P> {
    type Item = &'a P;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|&n| unsafe { &*n }.p())
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}
impl<'a, P> ExactSizeIterator for NeighboursData<'a, P> {}

pub struct IncidentEdges<'a, P: 'a> {
    iter: Neighbours<'a, P>,
    node: &'a Node<P>
}
impl<'a, P> Iterator for IncidentEdges<'a, P> {
    type Item = Edge<'a, P>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|n| Edge { n1: &self.node, n2: n})
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}
impl<'a, P> ExactSizeIterator for IncidentEdges<'a, P> {}

impl<P: Point> Node<P> {
    pub fn init(&mut self) {
        let mut p = self.p.clone();
        p.mul(P::R::zero());
        let num = P::R::from(self.neighbours().len() as f64);
        for n in self.neighbours() {
            p.add(n.p().clone().div(num));
        }
        self.p = p;
    }
}

impl<P> Node<P> {
    pub fn p(&self) -> &P {
        &self.p
    }

    pub fn p_mut(&mut self) -> &mut P {
        &mut self.p
    }

    pub fn id(&self) -> usize {
        self.id
    }

    pub fn is_terminal(&self) -> bool {
        self.is_terminal
    }

    pub fn edges(&self) -> IncidentEdges<P> {
        IncidentEdges {
            iter: self.neighbours(),
            node: &self
        }
    }

    pub fn neighbours(&self) -> Neighbours<P> {
        Neighbours {
            iter: self.ns.iter()
        }
    }

    pub fn neighbours_data_mut(&mut self) -> (&mut P, NeighboursData<P>) {
        (&mut self.p, NeighboursData {
            iter: self.ns.iter()
        })
    }

    pub fn neighbours_mut(&mut self) -> NeighboursMut<P> {
        NeighboursMut {
            iter: self.ns.iter()
        }
    }
}

impl<P: fmt::Display> fmt::Display for Node<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "{} [ ", self.p));
        try!(self.neighbours().fold(Ok(()), |_: fmt::Result, n|
            write!(f, "#{} ", n.id())));
        write!(f, "]")
    }
}
