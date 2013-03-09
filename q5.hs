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

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Axis.Floating
import Data.Accessor
import Data.Colour
import Data.Colour.Names
import Numeric

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

-- | Line search by Newton's method
newtonSearch :: (Num a) => LineSearch f a
newtonSearch f df p x = undefined

-- | Line search by secant method with given tolerance
secantSearch :: (Num a, Fractional a) => a -> LineSearch f a
secantSearch eps f df p x = undefined

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

-- | 'secant f x0 x1' is the series of secant method approximations of
-- a root of 'f' surrounded by starting points 'x0' and 'x1'
secant :: (Num a, Fractional a) => (a -> a) -> a -> a -> [a]
secant f = curry go
  where go (x0, x1) = let x2 = x1 - f x1 * (x1 - x0) / (f x1 - f x0)
                      in x2 : go (x1, x2)

-- | Newton's method
newton :: (Num a, Ord a, Additive f, Metric f)
       => (f a -> a) -> (f a -> f a) -> (f a -> f (f a)) -> f a -> [f a]
newton f df ddfInv x0 = iterate go x0
  where go x = x ^-^ ddfInv x !* df x

-- | Barzilai-Borwein 1988 is a non-monotonic optimization method
barzilaiBorwein :: (Additive f, Metric f, Functor f, Fractional a)
                => (f a -> a) -> (f a -> f a) -> f a -> f a -> [f a]
barzilaiBorwein f df x0 x1 = go (x0, x1)
  where go (x0,x1) = let s = x1 ^-^ x0
                         z = df x1 ^-^ df x0
                         alpha = (s `dot` z) / (z `dot` z)
                         x2 = x1 ^-^ alpha *^ df x1
                     in x2 : go (x1, x2)

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
     forM_ (take 10 $ conjGrad search fletcherReeves f df x0) $ \x->do print (x, f x)

     putStrLn "\n\nSteepest descent"
     forM_ (take 10 $ gradientDescent search f df x0) $ \x->do print (x, f x)

     putStrLn "\n\nNewton"
     let ddfInv = maybe (error "Can't invert Hessian") id
                  . inv22 . hessian rosenbrock
     forM_ (take 10 $ newton f df ddfInv x0) $ \x->do print (x, f x)

     let x0s = [ V2 0.01 2.8, V2 0 3, V2 0.9 0.9 ]
         logAxis = autoScaledLogAxis $ loga_labelf ^= (\x->showFFloat (Just 2) (log x) "") $ defaultLogAxis
         layout = layout1_plots ^= (map (Left . toPlot) $ plot search fletcherReeves x0s)
                  $ layout1_left_axis ^: laxis_generate ^= logAxis
                  $ defaultLayout1
     renderableToPDFFile (toRenderable layout) 600 600 "q5.pdf"

plotValue :: [V2 Double] -> PlotLines Int Double
plotValue xs =
    plot_lines_values ^= [take 50 $ zip [1..] $ map rosenbrock xs]
    $ defaultPlotLines

plotError :: V2 Double -> [V2 Double] -> PlotLines Int Double
plotError x xs =
    plot_lines_values ^= [take 50 $ zip [1..] $ map (x `qd`) xs]
    $ defaultPlotLines

plot :: LineSearch V2 Double -> Beta V2 Double -> [V2 Double] -> [PlotLines Int Double]
plot search beta x0s =
    let f = rosenbrock
        df = grad rosenbrock :: V2 Double -> V2 Double
        ddfInv = maybe (error "Can't invert Hessian") id
                 . inv22 . hessian rosenbrock
        showV2 (V2 x y) = "("++show x++","++show y++")"
        cg = map (\x0->("conj. grad. "++showV2 x0, blue, conjGrad search beta f df x0)) x0s
        gd = map (\x0->("grad. desc. "++showV2 x0, green, gradientDescent search f df x0)) x0s
        n  = map (\x0->("newton "++showV2 x0, red, newton f df ddfInv x0)) x0s
        errorPlots = map (\(title,color,xs) ->
                             plot_lines_title ^= title
                             $ plot_lines_style ^: line_color ^= opaque color
                             $ plotError (V2 1 1) $ take 50 xs
                         ) $ concat [cg, gd, n]
    in errorPlots
