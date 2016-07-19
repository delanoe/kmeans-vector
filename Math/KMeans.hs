{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE CPP #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE BangPatterns #-}

{- |
Module      :  Math.KMeans
Copyright   :  (c) Alp Mestanogullari, Ville Tirronen, 2011-2015
License     :  BSD3
Maintainer  :  Alp Mestanogullari <alpmestan@gmail.com>
Stability   :  experimental

An implementation of the k-means clustering algorithm based on the vector package.

The core functions of this module are 'kmeans' and 'kmeansWith'. See some examples
on <http://github.com/alpmestan/kmeans-vector github>.

-}
module Math.KMeans
  ( -- * The meat of this package: 'kmeans' 
    kmeans
  , kmeans'
  , kmeansWith
  , kmeansWith'

  , -- * Types
    Distance
  , Clusters
  , Cluster(..)
  , Centroids

  , -- * Misc.
    partition
  , euclidSq
  , l1dist
  , linfdist
  ) where

import Control.Monad.Identity
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.List as L
import Data.Function (on)

#if !MIN_VERSION_base(4, 8, 0)
import Data.Foldable
import Data.Monoid
import Data.Traversable
#endif

-- | A distance on vectors
type Distance = U.Vector Double -> U.Vector Double -> Double

-- | The euclidean distance without taking the final square root
--   This would waste cycles without changing the behavior of the algorithm
euclidSq :: Distance
euclidSq v1 v2 = U.sum $ U.zipWith diffsq v1 v2
  where diffsq a b = (a-b)^(2::Int)
{-# INLINE euclidSq #-}

-- | L1 distance of two vectors: d(v1, v2) = sum on i of |v1_i - v2_i|
l1dist :: Distance
l1dist v1 v2 = U.sum $ U.zipWith diffabs v1 v2
  where diffabs a b = abs (a - b)
{-# INLINE l1dist #-}

-- | L-inf distance of two vectors: d(v1, v2) = max |v1_i - v2_i]
linfdist :: Distance
linfdist v1 v2 = U.maximum $ U.zipWith diffabs v1 v2
  where diffabs a b = abs (a - b)
{-# INLINE linfdist #-}

-- | This is what 'kmeans' hands you back. It's just a 'V.Vector' of clusters
--   that will hopefully be of length 'k'.
type Clusters a = V.Vector (Cluster a)

-- | This type is used internally by 'kmeans'. It represents our (hopefully)
--   @k@ centroids, obtained by computing the new centroids of a 'Cluster'
type Centroids  = V.Vector (U.Vector Double)

-- | A 'Cluster' of points is just a list of points
newtype Cluster a =
  Cluster { elements :: [a] -- ^ elements that belong to that cluster
          } deriving (Eq, Show, Functor, Monoid, Foldable, Traversable)

clusterAdd :: Cluster a -> a -> Cluster a
clusterAdd (Cluster c) x = Cluster (x:c)

emptyCluster :: Cluster a
emptyCluster = Cluster []

addCentroids :: U.Vector Double -> U.Vector Double -> U.Vector Double
addCentroids v1 v2 = U.zipWith (+) v1 v2

-- | This is the current partitionning strategy used
--   by 'kmeans'. If we want @k@ clusters, we just 
--   try to regroup consecutive elements in @k@ buckets
partition :: Int -> [a] -> Clusters a
partition k vs = V.fromList $ go vs
  where go l = case L.splitAt n l of
          (vs', []) -> [Cluster vs']
          (vs', vss) -> Cluster vs' : go vss
        n = (length vs + k - 1) `div` k

-- | Run the kmeans clustering algorithm.
-- 
--  > kmeans f distance k points
-- 
-- will run the algorithm using 'f' to extract features from your type,
-- using 'distance' to measure the distance between vectors,
-- trying to separate 'points' in 'k' clusters.
--
-- Extracting features just means getting a 'U.Vector'
-- with 'Double' coordinates that will represent your type
-- in the space in which 'kmeans' will run.
kmeans :: (a -> U.Vector Double) -- ^ feature extraction
       -> Distance               -- ^ distance function
       -> Int                    -- ^ the 'k' to run 'k'-means with (i.e number of desired clusters)
       -> [a]                    -- ^ input list of 'points'
       -> Clusters a             -- ^ result, hopefully 'k' clusters of points
kmeans extract dist k points = 
  runIdentity $ kmeansWith (\n ps -> pure $ partition n ps) extract dist k points


kmeans' :: (a -> U.Vector Double) -- ^ feature extraction
       -> Distance               -- ^ distance function
       -> Int                    -- ^ the 'k' to run 'k'-means with (i.e number of desired clusters)
       -> [a]                    -- ^ input list of 'points'
       -> (Centroids, Clusters a)             -- ^ result, hopefully 'k' clusters of points
kmeans' extract dist k points = 
  runIdentity $ kmeansWith' (\n ps -> pure $ partition n ps) extract dist k points


-- | Same as 'kmeans', except that instead of using 'partition', you supply your own
--   function for choosing the initial clusterinV. Two important things to note:
-- 
--   * If you don't need any kind of effect and just have a 'partition'-like function
--     you want to use, @m@ will can just be 'Identity' here. If that's too 
--     obnoxious to work with, please let me know and I may just provide a separate
--     'kmeansWith' function with no there. But most of the time, you'll probably just
--     be interested in the following scenario.
-- 
--   * Most likely, you want to have something smarter than our simple 'partition' function.
--     A couple of papers I have read claim very decent results by using some precise
--     probabilistic schemas for the initial partitionninV. In this case, your @m@ would
--     probably be 'IO' or 'ST' (e.g using my <http://hackage.haskell.org/package/probable probable> package)
--     and you could fine-tune the way the initial clusters are picked so that the algorithm
--     may give better results. Of course, if your initialization is monadic, so is the result. 
kmeansWith :: forall m a. Monad m
           => (Int -> [a] -> m (Clusters a)) -- ^ how should we partition the points?
           -> (a -> U.Vector Double)         -- ^ get the coordinates of a "point"
           -> Distance                       -- ^ what distance do we use
           -> Int                            -- ^ number of desired clusters
           -> [a]                            -- ^ list of points
           -> m (Clusters a)                 -- ^ resulting clustering
kmeansWith initF extract dist k points = go `liftM` initF k points
  
  where 
    go :: Clusters a -> Clusters a
    go pgroups =
      case kmeansStep pgroups of
        pgroups' | pgroupsEqualUnder pgroups pgroups'  -> pgroups
                 | otherwise -> go pgroups' 

    kmeansStep :: Clusters a -> Clusters a
    kmeansStep clusters = 
      case centroidsOf clusters of
        centroids -> 
            V.filter (not . null . elements)
          . V.unsafeAccum clusterAdd (V.replicate k emptyCluster)
          . map (pairToClosestCentroid centroids)
          $ points

    centroidsOf :: Clusters a -> Centroids
    centroidsOf cs = V.map centroidOf cs
      where

        centroidOf (Cluster elts) = 
            U.map (/n) 
          . L.foldl1' addCentroids
          $ map extract elts

          where n = fromIntegral (length elts)

    pairToClosestCentroid :: Centroids -> a -> (Int, a)
    pairToClosestCentroid cs a = (minDistIndex, a)
      where !minDistIndex = V.minIndexBy (compare `on` dist (extract a)) cs

    pgroupsEqualUnder :: Clusters a -> Clusters a -> Bool
    pgroupsEqualUnder g1 g2 = elementsOf g1 == elementsOf g2
      where elementsOf = V.map (map extract . elements)

kmeansWith' :: forall m a. Monad m
           => (Int -> [a] -> m (Clusters a)) -- ^ how should we partition the points?
           -> (a -> U.Vector Double)         -- ^ get the coordinates of a "point"
           -> Distance                       -- ^ what distance do we use
           -> Int                            -- ^ number of desired clusters
           -> [a]                            -- ^ list of points
           -> m (Centroids, Clusters a)      -- ^ resulting clustering
kmeansWith' initF extract dist k points = go `liftM` initF k points

  where 
    go :: Clusters a -> (Centroids, Clusters a)
    go pgroups =
      case kmeansStep pgroups of
        (c, pgroups') | pgroupsEqualUnder pgroups pgroups'  -> (c, pgroups)
                      | otherwise -> go pgroups'

    kmeansStep :: Clusters a -> (Centroids, Clusters a)
    kmeansStep clusters =
      case centroidsOf clusters of
        centroids -> (centroids, clusters) where
            clusters = V.filter (not . null . elements)
                       . V.unsafeAccum clusterAdd (V.replicate k emptyCluster)
                       . map (pairToClosestCentroid centroids)
                       $ points

    centroidsOf :: Clusters a -> Centroids
    centroidsOf cs = V.map centroidOf cs
      where

        centroidOf (Cluster elts) =
            U.map (/n) 
          . L.foldl1' addCentroids
          $ map extract elts

          where n = fromIntegral (length elts)

    pairToClosestCentroid :: Centroids -> a -> (Int, a)
    pairToClosestCentroid cs a = (minDistIndex, a)
      where !minDistIndex = V.minIndexBy (compare `on` dist (extract a)) cs

    pgroupsEqualUnder :: Clusters a -> Clusters a -> Bool
    pgroupsEqualUnder g1 g2 = elementsOf g1 == elementsOf g2
      where elementsOf = V.map (map extract . elements)


{-# INLINE kmeansWith #-}
