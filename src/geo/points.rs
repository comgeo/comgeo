use std::fmt;
use std::slice::{IterMut, Iter};

use traits::*;

#[derive(Debug)]
pub struct PointNd<N> {
    coords: Vec<N>,
    id: usize
}

impl<N> PointNd<N> {
    pub fn new(c: Vec<N>) -> PointNd<N> {
        PointNd { coords: c, id: 0 }
    }
}

impl<R: Real> Point for PointNd<R> {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords[..].iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords[..].iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Real + Clone> Clone for PointNd<E> {
    fn clone(&self) -> Self {
        PointNd { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for PointNd<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: Real + fmt::Display> fmt::Display for PointNd<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "( "));
        try!(self.iter().fold(Ok(()), |_, n: &E| write!(f, "{} ", n)));
        write!(f, ")")
    }
}

#[derive(Debug)]
pub struct Point2d<E> {
    coords: [E; 2],
    id: usize
}

impl<E> Point2d<E> {
    pub fn new(c: [E; 2]) -> Point2d<E> {
        Point2d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 2] {
        &self.coords
    }
}

impl<R: Real> Point for Point2d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Real + Clone> PartialEq for Point2d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: Copy> Clone for Point2d<E> {
    fn clone(&self) -> Self {
        Point2d { id: self.id, coords: [self.coords[0], self.coords[1]] }
    }
}

impl<E: fmt::Display> fmt::Display for Point2d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({},{})", self.coords[0], self.coords[1])
    }
}

#[derive(Debug)]
pub struct Point3d<E> {
    coords: [E; 3],
    id: usize
}
#[derive(Debug)]
pub struct Point4d<E> {
    coords: [E; 4],
    id: usize
}
#[derive(Debug)]
pub struct Point5d<E> {
    coords: [E; 5],
    id: usize
}
#[derive(Debug)]
pub struct Point6d<E> {
    coords: [E; 6],
    id: usize
}
#[derive(Debug)]
pub struct Point7d<E> {
    coords: [E; 7],
    id: usize
}
#[derive(Debug)]
pub struct Point8d<E> {
    coords: [E; 8],
    id: usize
}
#[derive(Debug)]
pub struct Point9d<E> {
    coords: [E; 9],
    id: usize
}

impl<E> Point3d<E> {
    pub fn new(c: [E; 3]) -> Point3d<E> {
        Point3d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 3] {
        &self.coords
    }
}

impl<R: Real> Point for Point3d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point3d<E> {
    fn clone(&self) -> Self {
        Point3d { id: self.id, coords: [self.coords[0], self.coords[1], self.coords[2]] }
    }
}

impl<E: Real + Clone> PartialEq for Point3d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point3d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({},{},{})", self.coords[0], self.coords[1], self.coords[2])
    }
}


impl<E> Point4d<E> {
    pub fn new(c: [E; 4]) -> Point4d<E> {
        Point4d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 4] {
        &self.coords
    }
}

impl<R: Real> Point for Point4d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point4d<E> {
    fn clone(&self) -> Self {
        Point4d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point4d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point4d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}


impl<E> Point5d<E> {
    pub fn new(c: [E; 5]) -> Point5d<E> {
        Point5d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 5] {
        &self.coords
    }
}

impl<R: Real> Point for Point5d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point5d<E> {
    fn clone(&self) -> Self {
        Point5d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point5d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point5d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}

impl<E> Point6d<E> {
    pub fn new(c: [E; 6]) -> Point6d<E> {
        Point6d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 6] {
        &self.coords
    }
}

impl<R: Real> Point for Point6d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point6d<E> {
    fn clone(&self) -> Self {
        Point6d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point6d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point6d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}

impl<E> Point7d<E> {
    pub fn new(c: [E; 7]) -> Point7d<E> {
        Point7d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 7] {
        &self.coords
    }
}

impl<R: Real> Point for Point7d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point7d<E> {
    fn clone(&self) -> Self {
        Point7d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point7d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point7d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}

impl<E> Point8d<E> {
    pub fn new(c: [E; 8]) -> Point8d<E> {
        Point8d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 8] {
        &self.coords
    }
}

impl<R: Real> Point for Point8d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point8d<E> {
    fn clone(&self) -> Self {
        Point8d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point8d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point8d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}

impl<E> Point9d<E> {
    pub fn new(c: [E; 9]) -> Point9d<E> {
        Point9d { coords: c, id: 0 }
    }

    pub fn arr(&self) -> &[E; 9] {
        &self.coords
    }
}

impl<R: Real> Point for Point9d<R>  {
    type R = R;

    fn id(&self) -> usize {
        self.id
    }

    fn set_id(&mut self, id: usize) {
        self.id = id;
    }

    fn dim(&self) -> usize {
        self.coords.len()
    }

    fn iter(&self) -> Iter<R> {
        self.coords.iter()
    }

    fn iter_mut(&mut self) -> IterMut<R> {
        self.coords.iter_mut()
    }

    fn coords(&self) -> &[R] {
        &self.coords[..]
    }

    fn coords_mut(&mut self) -> &mut [R] {
        &mut self.coords[..]
    }
}

impl<E: Copy> Clone for Point9d<E> {
    fn clone(&self) -> Self {
        Point9d { id: self.id, coords: self.coords.clone() }
    }
}

impl<E: Real + Clone> PartialEq for Point9d<E> {
    fn eq(&self, other: &Self) -> bool {
        self.iter().zip(other.iter()).all(&|(c, oc)| c == oc)
    }
}

impl<E: fmt::Display> fmt::Display for Point9d<E> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "({}", self.coords[0]));
        for c in self.coords[1..].iter() {
            try!(write!(f, ",{}", c));
        }
        write!(f, ")")
    }
}
