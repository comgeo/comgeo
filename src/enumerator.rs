use traits::*;
use steinertree::*;
use algorithms::mst::*;

use std::marker::PhantomData;
use std::fmt;
use std::iter::*;
use std::time::{Duration, Instant};
use std::ops::{Add};
use std::io::{self, BufWriter, Write};

#[derive(PartialEq, Debug)]
enum GPState {
    Done,
    Start,
    Running
}

pub struct GPEnumeration<P: Point, S: TerminalSorter> {
    tree: SteinerTree<P>,
    edges: Vec<(usize, usize)>,
    top: Vec<usize>,
    t: Vec<P>,
    bsd: Vec<Vec<P::R>>,
    ss: Vec<P::R>,
    t_len: usize,
    sorter: S,

    e_bsd: bool,
    e_ss: bool,
    state: GPState,
    data: GPEnumerationData,
    _m: PhantomData<P>
}

impl<P: Point, S: TerminalSorter> fmt::Debug for GPEnumeration<P, S> {
    fn fmt(&self, _: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        Ok(())
    }
}

impl<P: Point> Default for GPEnumeration<P, FurthestSiteOrdering> {
    fn default() -> Self {
        GPEnumeration::new(true, true, FurthestSiteOrdering)
    }
}


pub struct NoOrdering;
impl TerminalSorter for NoOrdering {
    fn sort<P, M>(&mut self, _: &mut[P], _: &M)
        where P: Point, M: MinkowskiSpace<P> { }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "no ordering")
    }
}
impl fmt::Display for NoOrdering {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}


pub struct FurthestSiteOrdering;
impl TerminalSorter for FurthestSiteOrdering {
    fn sort<P, M>(&mut self, t: &mut[P], geo: &M)
        where P: Point, M: MinkowskiSpace<P> {

        let mut sorting = [0; 3];

        // Three with max sum of distances
        let mut max = P::R::zero();
        for i in 0..t.len()-2 {
            for j in i+1..t.len()-1 {
                let dij = geo.dist(&t[i], &t[j]);
                for k in j+1..t.len() {
                    let dik = geo.dist(&t[i], &t[k]);
                    let dkj = geo.dist(&t[k], &t[j]);
                    let sum = dij + dik + dkj;
                    if sum > max {
                        max = sum;
                        sorting[0] = i;
                        sorting[1] = j;
                        sorting[2] = k;
                    }
                }
            }
        }
        t.swap(0, sorting[0]);
        t.swap(1, sorting[1]);
        t.swap(2, sorting[2]);

        // Add rest by max distance to any previous added
        for i in 3..t.len()-1 {
            let k = {
                let (sorted, rest) = t.split_at(i);
                rest.iter().enumerate().fold((0, P::R::zero()),
                    |(ai, max), (i, p)| {
                        let maxdist = sorted.iter().fold(P::R::zero(),
                            |max, a| max.max(geo.dist(p, a)));
                        if maxdist > max {
                            (i, maxdist)
                        } else {
                            (ai, max)
                        }
                    }
                ).0
            };
            t.swap(i, i+k);
        }
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "furthest site ordering")
    }
}
impl fmt::Display for FurthestSiteOrdering {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}

impl<P: Point, S: TerminalSorter> GPEnumeration<P, S> {
    fn new(e_bsd: bool, e_ss: bool, sorter: S) -> Self {
        GPEnumeration {
            e_bsd: e_bsd,
            e_ss: e_ss,
            state: GPState::Done,
            tree: SteinerTree::new(&[], &[], &[]),
            edges: Vec::new(),
            top: Vec::new(),
            t: Vec::new(),
            sorter: sorter,
            bsd: Vec::new(),
            ss: Vec::new(),
            t_len: 0,
            data: GPEnumerationData::new(0),
            _m: PhantomData
        }
    }

    pub fn default_with_sorter(sorter: S) -> Self {
        GPEnumeration::new(true, true, sorter)
    }

    #[inline]
    fn get_bottleneck(&self, p1: &P, p2: &P) -> P::R {
        self.bsd[p1.id()][p2.id()]
    }

    #[inline]
    fn get_ss(&self, p: &P) -> P::R {
        self.ss[p.id()]
    }

    fn calc_bsd<M>(&mut self, terms: &[P], geo: &M)
        where M: MinkowskiSpace<P> {

        fn dfs_bottleneck<M, P>(start: &Node<P>, goal: &Node<P>, prev: &Node<P>, b: P::R, geo: &M) -> (bool, P::R)
            where P: Point, M: MinkowskiSpace<P> {

            for e in start.edges().filter(|e| e.nodes().1.p() != prev.p()) {
                if e.nodes().1.p() == goal.p() {
                    return (true, b.max(e.len(geo)));
                }
                if let res @ (true, _) = dfs_bottleneck(&e.nodes().1, goal, start, b.max(e.len(geo)), geo) {
                    return res;
                }
            }
            (false, P::R::zero())
        }

        let mst = Kruskal::new().find(terms, geo);
        self.bsd = vec![vec![P::R::zero(); terms.len()]; terms.len()];

        for (n1, n2) in mst.node_pairs() {
            let (_, b) = dfs_bottleneck(&n1, &n2, &n1, n1.edges().next().unwrap().len(geo), geo);
            self.bsd[n1.p().id()][n2.p().id()] = b;
            self.bsd[n2.p().id()][n1.p().id()] = b;
        }
    }

    fn calc_ss<M>(&mut self, terms: &[P], geo: &M)
        where M: MinkowskiSpace<P> {

        self.ss = vec![P::R::zero(); terms.len()];
        for t in terms {
            let mut iter = terms.iter().filter(|p| p.id() != t.id());
            let mut d = geo.dist(iter.next().unwrap(), t);
            for n in iter {
                let dd = geo.dist(n, t);
                if dd < d {
                    d = dd;
                }
            }
            self.ss[t.id()] = d;
        }
    }

    fn pop(&mut self) {
        self.top.pop();
        let ei = *self.top.last().unwrap() - 1;
        let (b, si) = self.edges.pop().unwrap();
        let (ti, _) = self.edges.pop().unwrap();
        let a = self.edges[ei].0;
        self.edges[ei] = (a, b);
        self.tree.add_edge((a, b));

        self.tree.pop_node();
        let tp = self.tree.pop_node().unwrap();
        self.t.push(tp);
    }

    fn backtrack(&mut self) -> bool {
        if self.t.len() == 0 {
            self.pop();
        }

        while *self.top.last().unwrap() == self.edges.len()  {
            if self.top.len() == 1 {
                self.state = GPState::Done;
                return false;
            } else {
                self.pop();
            }
        }

        true
    }

    fn insert_on_edge(&mut self, ei: usize) -> (usize, usize) {
        //println!("Insert_on_edge: start");
        //stdout().flush();
        let p = self.t.pop().unwrap();
        let (a, b) = self.edges[ei];
        self.tree.remove_edge((a, b));
        let steiner = p.clone();
        let ti = self.tree.push_node(p, &[], true).id();
        let si = self.tree.push_node(steiner, &[a, b, ti], false).id();

        self.edges[ei] = (a, si);
        self.edges.push((ti, si));
        self.edges.push((b, si));

        (ti, si)
    }


    fn prune_paths_check<F>(s: &Node<P>, t: &Node<P>, f: &mut F) -> bool
        where F: FnMut(&P, &P, usize) -> bool {

        fn aux<F, P>(cur: &Node<P>, prev: &Node<P>, p: &P, k: usize, f: &mut F) -> bool
            where F: FnMut(&P, &P, usize) -> bool {

            if cur.is_terminal() {
                f(cur.p(), p, k)
            } else {
                cur.neighbours()
                   .filter(|n| n.id() != prev.id())
                   .any(|n| aux(n, cur, p, k + 1, f))
            }
        }

        let mut iter = s.neighbours().filter(|n| n.id() != t.id());
        let n1 = iter.next().unwrap();
        let n2 = iter.next().unwrap();
        aux(n1, s, t.p(), 1, &mut |ps, p, k| f(ps, p, k) || (aux(n2, s, ps, k, f)))
            || aux(n2, s, t.p(), 1, f)
    }

    fn prune_check<M>(&mut self, ti: usize, si: usize, geo: &M) -> bool
        where M: MinkowskiSpace<P> {

        let start = Instant::now();
        let depth = self.t_len - self.t.len() - 4;

        let prune = {
            // Safe if we do NOT use self.data afterwards.
            let data = unsafe { &mut *(&mut self.data as *mut GPEnumerationData) };
            let (t, s) = (self.tree.i(ti), self.tree.i(si));
            if self.e_bsd || self.e_ss {
                Self::prune_paths_check(s, t, &mut |p1, p2, k| {
                    let j = self.t.len();
                    let d = geo.dist(p1, p2);
                    let b = self.get_bottleneck(p1, p2);
                    if self.e_bsd && d > P::R::from((k+j+1) as f64) * b {
                        data.bsd_pruned[depth] += 1;
                        true
                    } else if self.e_ss && d > self.get_ss(p1) + self.get_ss(p2)
                                               + P::R::from((k+j-1) as f64) * b
                    {
                        data.ss_pruned[depth] += 1;
                        true
                    } else {
                        false
                    }
                })
            } else {
                false
            }
        };

        self.data.prune_times[depth] += Instant::now() - start;
        prune
    }
}

#[derive(Debug, Clone)]
pub struct GPEnumerationData {
    nodes: usize,
    time: Duration,
    bsd_pruned: Vec<usize>,
    ss_pruned: Vec<usize>,
    prune_times: Vec<Duration>,
    sort_time: Duration,
    bsd_time: Duration,
    ss_time: Duration,
    init_time: Duration,
}

impl GPEnumerationData {
    fn new(mut ts: usize) -> Self {
        ts = if ts < 3 { 0 } else { ts - 3 };
        GPEnumerationData {
            nodes: 0,
            time: Duration::new(0, 0),
            bsd_pruned: vec![0; ts],
            ss_pruned: vec![0; ts],
            prune_times: vec![Duration::new(0, 0); ts],
            sort_time: Duration::new(0, 0),
            bsd_time: Duration::new(0, 0),
            ss_time: Duration::new(0, 0),
            init_time: Duration::new(0, 0)
        }
    }

    fn clear(&mut self, mut ts: usize) {
        ts = if ts < 3 { 0 } else { ts - 3 };
        self.nodes = 0;
        self.time = Duration::new(0,0);
        self.bsd_pruned = vec![0; ts];
        self.ss_pruned = vec![0; ts];
        self.prune_times = vec![Duration::new(0, 0); ts];
        self.sort_time = Duration::new(0,0);
        self.bsd_time = Duration::new(0,0);
        self.ss_time = Duration::new(0,0);
        self.init_time = Duration::new(0,0);
    }

    pub fn ss_pruned(&self) -> &[usize]  {
        &self.ss_pruned
    }

    pub fn bsd_pruned(&self) -> &[usize]  {
        &self.bsd_pruned
    }

    pub fn prune_times(&self) -> &[Duration]  {
        &self.prune_times
    }

    pub fn sort_time(&self) -> &Duration  {
        &self.sort_time
    }

    pub fn bsd_time(&self) -> &Duration  {
        &self.bsd_time
    }

    pub fn ss_time(&self) -> &Duration  {
        &self.ss_time
    }

    pub fn init_time(&self) -> &Duration  {
        &self.init_time
    }
}

impl EnumeratorData for GPEnumerationData {
    fn time(&self) -> &Duration {
        &self.time
    }

    fn nodes(&self) -> usize  {
        self.nodes
    }

    fn pruned(&self) -> usize  {
        self.bsd_pruned.iter().fold(0, Add::add) +
        self.ss_pruned.iter().fold(0, Add::add)
    }
}

impl fmt::Display for GPEnumerationData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Data for the Gilbert-Pollak enumeration:"));
        try!(writeln!(f, "\tNumber of nodes enumerated (i.e. not pruned): {}", self.nodes));

        try!(writeln!(f, "\n\tInit times:"));
        try!(writeln!(f, "\t\tBottleneck Steiner distances: {}", printdur(&self.bsd_time)));
        try!(writeln!(f, "\t\tSmallest sphere distances: {}", printdur(&self.ss_time)));
        try!(writeln!(f, "\t\tSorting: {}", printdur(&self.sort_time)));
        try!(writeln!(f, "\t\tTotal init time: {}", printdur(&self.init_time)));

        try!(writeln!(f, "\n\tNumber of pruned nodes & time used on pruning at each level:"));
        try!(self.prune_times.iter().zip(self.bsd_pruned.iter().zip(self.ss_pruned.iter())).enumerate()
            .fold(Ok(()), |_: fmt::Result, (level, (time, (bsd, ss)))|
                writeln!(f, "\t\tLevel {} - Time: {}, bsd pruned: {}, ss pruned: {}",
                    level+1, printdur(time), bsd, ss)));
        Ok(())
    }
}

impl<P: Point, S: TerminalSorter> Enumerator<P> for GPEnumeration<P, S> {
    type D = GPEnumerationData;

    fn init<M>(&mut self, mut terms: Vec<P>, geo: &M)
        where M: MinkowskiSpace<P> {

        self.data.clear(terms.len());

        let start = Instant::now();
        self.sorter.sort(&mut terms[..], geo);
        self.data.sort_time = Instant::now() - start;

        self.top.clear();
        self.t.clear();
        self.edges.clear();
        self.t_len = terms.len();
        self.state = GPState::Start;

        for (i, t) in terms.iter_mut().enumerate() {
            t.set_id(i);
        }

        if self.e_bsd {
            let bsdstart = Instant::now();
            self.calc_bsd(&terms[..], geo);
            self.data.bsd_time = Instant::now() - bsdstart;
        }

        if self.e_ss {
            let ssstart = Instant::now();
            self.calc_ss(&terms[..], geo);
            self.data.ss_time = Instant::now() - ssstart;
        }

        if terms.len() < 3 {
            self.tree = SteinerTree::new(&terms[..], &[], &[]);
            if terms.len() == 2 {
                self.tree.add_edge((0,1));
                self.edges.push((0,1));
            }
        } else {
            let p = terms[0].clone();
            self.tree = SteinerTree::new(&terms[self.t_len-3..], &[p], &[(0,3), (1,3), (2,3)]);
            self.tree.steiner_i(3).init();
            terms.truncate(self.t_len - 3);

            self.edges.extend_from_slice(&[(0,3),(1,3),(2,3)]);
            self.top.push(0);
            self.t = terms;
        }
        self.data.init_time = Instant::now() - start;
    }

    fn next<M>(&mut self, geo: &M) -> bool
        where M: MinkowskiSpace<P> {

        let start = Instant::now();

        if self.state == GPState::Start {
            if self.t.len() == 0 {
                self.state = GPState::Done;
            } else {
                self.state = GPState::Running;
            }
            self.data.time = Instant::now() - start;
            return true;
        }

        if self.state == GPState::Done || !self.backtrack() {
            self.data.time = Instant::now() - start;
            return false;
        }

        loop {
            let edge = *self.top.last().unwrap();
            let (ti, si) = self.insert_on_edge(edge);
            *self.top.last_mut().unwrap() += 1;
            self.top.push(0);

            //println!("Top: {:?}", self.top);
            //println!("Edges: {:?}", self.edges);
            //println!("Tree: {}", self.tree);
            if self.prune_check(ti, si, geo) {
                self.pop();
                if !self.backtrack() {
                    self.data.time = Instant::now() - start;
                    return false;
                }
            } else {
                self.data.nodes += 1;
                self.data.time = Instant::now() - start;
                return true;
            }
        }
    }



    fn tree(&self) -> &SteinerTree<P> {
        &self.tree
    }

    fn tree_mut(&mut self) -> &mut SteinerTree<P> {
        &mut self.tree
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        fn indent(f: &mut fmt::Formatter, indent: u32) -> fmt::Result {
            for _ in 0..indent {
                try!(write!(f, " "));
            }
            Ok(())
        }

        try!(writeln!(f, "Gilbert-Pollak enumerator (depth-first) with:"));
        if self.e_bsd {
            try!({indent(f, inde); writeln!(f, "  - bottleneck Steiner distance pruning;") });
        }
        if self.e_ss {
            try!({indent(f, inde); writeln!(f, "  - smallest spheres pruning;") });
        }

        try!({indent(f, inde); write!(f, "  - and "); self.sorter.print(f, inde+10) });
        Ok(())
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data())
    }
}

impl<P: Point, S: TerminalSorter> fmt::Display for GPEnumeration<P, S> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}
