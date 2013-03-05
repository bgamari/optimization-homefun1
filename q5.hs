import Data.Foldable (Foldable)
import qualified Data.Foldable as F
import Data.Traversable
import Data.Distributive (Distributive)
import Data.Functor.Bind (Apply)
import Control.Applicative
import Control.Monad
import Numeric.AD hiding (gradientDescent)
import Numeric.AD.Types (auto)
import qualified Numeric.LinearAlgebra as LA
import Linear

-- | A 'LineSearch' method 'search f df p x' determines a step size
-- in direction 'p' from point 'x' for function 'f' with gradient 'df'
type LineSearch f a = (f a -> a) -> (f a -> f a) -> f a -> f a -> a

-- | Armijo condition    
armijo :: (Num a, Additive f, Ord a, Metric f)
       => a -> (f a -> a) -> (f a -> f a) -> f a -> f a -> a -> Bool
armijo c1 f df x p a =
    f (x ^+^ a *^ p) <= f x + c1 * a * (df x `dot` p)

-- | Curvature condition
curvature :: (Num a, Ord a, Additive f, Metric f)
          => a -> (f a -> a) -> (f a -> f a) -> f a -> f a -> a -> Bool
curvature c2 f df x p a =
    df (x ^+^ a *^ p) `dot` p >= c2 * (df x `dot` p)

-- | Backtracking line search algorithm
backtrackingSearch :: (Num a, Ord a, Metric f)
                   => a -> a -> LineSearch f a
backtrackingSearch gamma c f df p x =
    head $ dropWhile (not . armijo c f df x p) $ iterate (*gamma) c

-- | Constant line search
constantSearch :: a -> LineSearch f a
constantSearch c f df p x = c

-- | Gradient descent method
gradientDescent :: (Num a, Ord a, Additive f, Metric f)
                => LineSearch f a -> (f a -> a) -> (f a -> f a) -> f a -> [f a]
gradientDescent search f df x0 = iterate go x0
  where go x = let p = negated (df x)
                   a = search f df p x
               in x ^+^ a *^ p
                
-- | A beta expression 'beta df0 df1 p' is an expression for the
-- conjugate direction contribution given the derivative 'df0' and
-- direction 'p' for iteration 'k', 'df1' for iteration 'k+1'
type Beta f a = f a -> f a -> f a -> a

-- | Conjugate gradient method with given beta and line search method
conjGrad :: (Num a, RealFloat a, Additive f, Metric f)
         => LineSearch f a -> Beta f a
         -> (f a -> a) -> (f a -> f a) -> f a -> [f a]
conjGrad search beta f df x0 = go (negated $ df x0) x0
  where go p x = let a = search f df p x
                     x' = x ^+^ a *^ p
                     b = beta (df x) (df x') p
                     p' = negated (df x') ^+^ b *^ p
                 in x' : go p' x'

-- | Fletcher-Reeves expression for beta
fletcherReeves :: (Num a, RealFloat a, Metric f) => Beta f a
fletcherReeves df0 df1 p0 = norm df1 / norm df0

-- | Polak-Ribiere expression for beta
polakRibiere :: (Num a, RealFloat a, Metric f) => Beta f a
polakRibiere df0 df1 p0 = df1 `dot` (df1 ^-^ df0) / norm df0

-- | Hestenes-Stiefel expression for beta
hestenesStiefel :: (Num a, RealFloat a, Metric f) => Beta f a
hestenesStiefel df0 df1 p0 =
    - (df1 `dot` (df1 ^-^ df0)) / (p0 `dot` (df1 ^-^ df0))

-- | Moore-Penrose pseudo-inverse
pseudoInverse :: (Functor m, Distributive n, Conjugate a)
              => m (n a) -> m (n a)
pseudoInverse a = undefined

-- | Inverse by iterative method of Ben-Israel and Cohen
-- with given starting condition
bicInv' :: (Functor m, Distributive m, Additive m,
            Applicative m, Apply m, Foldable m, Conjugate a)
        => m (m a) -> m (m a) -> [m (m a)]
bicInv' a0 a = iterate go a0
  where go ak = 2 *!! ak !-! ak !*! a !*! ak

-- | Inverse by iterative method of Ben-Israel and Cohen                
-- starting from 'alpha A^T'. Alpha should be set such that
-- 0 < alpha < 2/sigma^2 where sigma denotes the largest singular
-- value of A
bicInv :: (Functor m, Distributive m, Additive m,
           Applicative m, Apply m, Foldable m, Conjugate a)
       => a -> m (m a) -> [m (m a)]
bicInv alpha a = bicInv' (alpha *!! adjoint a) a
       
-- | Entry-wise matrix subtraction
(!-!) :: (Applicative m, Additive n, Num a)
      => m (n a) -> m (n a) -> m (n a)
a !-! b = (^-^) <$> a <*> b
infixl 6 !-!

-- | Newton's method
newton :: (Num a, Ord a, Additive f, Metric f)
       => (f a -> a) -> (f a -> f a) -> (f a -> f (f a)) -> f a -> [f a]
newton f df ddfInv x0 = iterate go x0
  where go x = x ^-^ ddfInv x !* df x
       
-- | Convert a nested-functor matrix to an hmatrix Matrix
toHMatrix :: (Functor m, Foldable m, Foldable n, LA.Element a)
          => m (n a) -> LA.Matrix a
toHMatrix = LA.fromLists . F.toList . fmap F.toList

-- | Rosenbrock function
rosenbrock :: Num a => V2 a -> a
rosenbrock (V2 x y) = (1-x)^2 + 100*(y-x^2)^2

main :: IO ()
main = do
     let f = rosenbrock
         df = grad rosenbrock :: V2 Double -> V2 Double
         x0 = V2 0.01 2.8
         --x0 = V2 0 3
         --x0 = V2 0.9 0.9
     let search = backtrackingSearch 0.1 0.2
     putStrLn "\n\nConjugate gradient"
     forM_ (take 100 $ conjGrad search fletcherReeves f df x0) $ \x->do print (x, f x)

     putStrLn "\n\nSteepest descent"
     forM_ (take 100 $ gradientDescent search f df x0) $ \x->do print (x, f x)

     putStrLn "\n\nNewton"
     --let ddfInv = head . take 100 . bicInv 0.1 . hessian rosenbrock
     let ddfInv = maybe (error "Can't invert Hessian") id
                  . inv22 . hessian rosenbrock
     forM_ (take 100 $ newton f df ddfInv x0) $ \x->do print (x, f x)
