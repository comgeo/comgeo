pub mod gilbert_pollak;

pub trait Enumerator<P: Point>: fmt::Display {
    type D: EnumeratorData;

    fn init<M: MinkowskiSpace<P>>(&mut self, Vec<P>, &M);
    fn next<M: MinkowskiSpace<P>>(&mut self, &M) -> bool;
    fn tree(&self) -> &SteinerTree<P>;
    fn tree_mut(&mut self) -> &mut SteinerTree<P>;
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}


pub trait EnumeratorData: fmt::Display + Clone {
    fn nodes(&self) -> usize;
    fn pruned(&self) -> usize;
    fn time(&self) -> &Duration;
}