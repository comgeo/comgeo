pub mod steiner_tree;
pub mod upper_bound;
pub mod steiner_bnb;
pub mod rmt;



pub trait SMT<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: SmtData;

    fn find(&mut self, Vec<P>, geo: &M) -> SteinerTree<P>;
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}


pub trait SmtData: fmt::Display + Clone {
    fn time(&self) -> &Duration;
    fn best_updates(&self) -> u64;
}