;;TAKEN FORM STACK OVERFLOW TO QUIET THOSE WARNINGS
(locally    (declare #+sbcl(sb-ext:muffle-conditions sb-kernel:redefinition-warning))
  (handler-bind
      (#+sbcl(sb-kernel:redefinition-warning #'muffle-warning))
    ;; stuff that emits redefinition-warning's
    ))

;;; Project 4: Build a simple evolutionary computation system.
(setf *random-state* (make-random-state t))
(require :sb-sprof)

;;; This is a nontrivial project to start with.  Do not discuss
;;; programming issues related to this project with other groups.
;;; If you have questions, contact me directly.  You may post to
;;; the Wiki questions of a general nature that you'd like people
;;; to respond to.

;;; Unquestionably the hardest code is the GP crossover operator and
;;; the PTC2 code.  For information about the Symbolic Regression
;;; and Artificial Ant problems, I have posted some chunks of my
;;; thesis to the web page to peruse.

;;; Please don't deviate much from this file -- I'm just looking for
;;; you to fill out these functions below and maybe add a bit more.  Rearranging
;;; things into multiple files etc. may be convenient for you but it's MUCH
;;; tougher to grade.

#|
In this project you will produce three things:

1. A high-level evolutionary computation framework

2. The representation, breeding functions, and evaluations to do a simple
GA for boolean vectors and for floating-point vectors, including a few test
functions on each.

3. The representation, breeding functions, and evaluations to do GP for
two problems:
A. Symbolic Regression
B. Artificial Ant

The high-level EC system will work as follows:

- Simple generational evolution
- GA-style Tournament selection
- A simple breeding function
- some simple statistics functions

I have provided the approximate function templates I myself used to complete
the task; with permission you can use these or go your own way, but otherwise
please try not to deviate much from this template.

The project is due approximately two and a half weeks from now or so.  Please 
try to get it in on time.

WHAT YOU MUST PROVIDE:

1. Completed code which works and compiles.  As simple as possible would be nice.

2. A short report describing the code you wrote and some experimentation with it.
Some things to try:
   -- different problems (in the vector section)
   -- different settings of mutation and crossover parameters
   -- different population sizes or run lengths
Try to get a handle on how problem difficulty changes when you tweak various
parameters.  Can you get mutation and crossover parameters which seem optimal for
a given problem for example?  (Obviously bigger population sizes are more optimal
but that's kinda cheating).  Note that this is a STOCHASTIC problem so you'll need
to run a number of times and get the mean best result in order to say anything of
consequence.  It'd be nice if the report were in LaTeX and it'd be good practice for
you as well, but it's not necessary.  I do not accept reports in Word.  Only send
me a PDF.


I VERY STRONGLY URGE YOU TO COMPILE YOUR CODE AND NOT RUN IT INTERPRETED!  If you're
running SBCL, it automatically compiles the code.  Else you'll need to compile it
manually.

Information on compiling code and doing compiler optimizations can be found in the
"Speed" chapter of Graham.
|#


;;; Useful Functions and Macros
(defparameter *debug* nil)

(defun dprint (some-variable &optional (additional-message '()))
	"Debug Print - useful for allowing error/status messages
to be printed while debug=t."
  (if *debug*
    (progn
      (if additional-message (print additional-message) nil)
      (print some-variable))
    some-variable))

(defmacro swap (elt1 elt2)
  "Swaps elt1 and elt2, using SETF.  Returns nil."
  (let ((temp (gensym)))
    `(let ((,temp ,elt1))
       (setf ,elt1 ,elt2)
       (setf ,elt2 ,temp)
       nil)))

(defmacro while (test return-value &rest body)
  "Repeatedly executes body as long as test returns true.
Then evaluates and returns return-value"
  `(progn (loop while ,test do (progn ,@body)) ,return-value))

(defun random-elt (sequence)
  "Returns a random element from sequence"
  (elt sequence (random (length sequence))))

(defun random? (&optional (prob 0.5))
  "Tosses a coin of prob probability of coming up heads,
then returns t if it's heads, else nil."
  (< (random 1.0) prob))

(defun generate-list (num function &optional no-duplicates)
  "Generates a list of size NUM, with each element created by
  (funcall FUNCTION).  If no-duplicates is t, then no duplicates
are permitted (FUNCTION is repeatedly called until a unique
new slot is created).  EQUALP is the test used for duplicates."
  (let (bag)
    (while (< (length bag) num) bag
      (let ((candidate (funcall function)))
	(unless (and no-duplicates
		     (member candidate bag :test #'equalp))
	  (push candidate bag))))))


;;;;;; TOP-LEVEL EVOLUTIONARY COMPUTATION FUNCTIONS 


;;; TOURNAMENT SELECTION

(defparameter *tournament-size* 7)

(defun tournament-select-one (population fitnesses)
  "Does one tournament selection and returns the selected individual. Algorithm 32."
	(let ((numbers ())
				(best-idx 0))
		(dprint fitnesses "TSO: ")
		(dotimes (p (1- *tournament-size*)) (setf numbers (append numbers (list (random (length population))))))
		(dprint numbers "Numbers: ")
		(setf best-idx (random (length population)))
		(dolist (next-idx numbers)
			(if (> (nth next-idx fitnesses) (nth best-idx fitnesses))
				(setf best-idx next-idx)))
		(dprint best-idx "Best-idx:")
		(nth best-idx population)
	)
)

(defun tournament-selector (num population fitnesses)
  "Does NUM tournament selections, and puts them all in a list"
	(let (chosen)
		;;(dprint num "Generations: ")
		;;(dprint population "Population: ")
		;;(dprint fitnesses "Fitnesses: ")
		(dotimes (w num)
			(setf chosen (append chosen (list (tournament-select-one population fitnesses)))))
		(dprint chosen "Chosen :")
	)
)

(defun simple-printer (pop fitnesses)  ;; I'm nice and am providing this for you.  :-)
  "Determines the individual in pop with the best (highest) fitness, then
prints that fitness and individual in a pleasing manner."
(let (best-ind best-fit)
    (print "random guy:")
    (print (get-random-element pop)) 
  
   (mapcar #'(lambda (ind fit)
		(when (or (not best-ind)
			  (< best-fit fit))
		  (setq best-ind ind)
		  (setq best-fit fit))) pop fitnesses)
    (format t "~%Best Individual of Generation...~%Fitness: ~a~%Individual:~a~%"
	    best-fit best-ind)
    fitnesses))

(defun evolve (generations pop-size
	       &key setup creator selector modifier evaluator printer)
  "Evolves for some number of GENERATIONS, creating a population of size
POP-SIZE, using various functions (Algorithm 20)."
  ;; The functions passed in are as follows:
  ;;(SETUP)                     called at the beginning of evolution, to set up
  ;;                            global variables as necessary
  ;;(CREATOR)                   creates a random individual
  ;;(SELECTOR num pop fitneses) given a population and a list of corresponding fitnesses,
  ;;                            selects and returns NUM individuals as a list.
  ;;                            An individual may appear more than once in the list.
  ;;(MODIFIER ind1 ind2)        modifies individuals ind1 and ind2 by crossing them
  ;;                            over and mutating them.  Returns the two children
  ;;                            as a list: (child1 child2).  Nondestructive to
  ;;                            ind1 and ind2.
  ;;(PRINTER pop fitnesses)     prints the best individual in the population, plus
  ;;                            its fitness, and any other interesting statistics
  ;;                            you think interesting for that generation.
  ;;(EVALUATOR individual)      evaluates an individual, and returns its fitness.
  ;: Pop will be guaranteed to be a multiple of 2 in size.
  ;;
  ;; HIGHER FITNESSES ARE BETTER

  ;; your function should call PRINTER each generation, and also print out or the
  ;; best individual discovered over the whole run at the end, plus its fitness
  ;; and any other statistics you think might be nifty.
;;; DAVID - could the variable declarations (below) go into the unused function "boolean-vector-sum-setup"? 
	(let ((population ())
				(fitnesses ()) ;; Might be able to use a lambda function to "create" fitnesses on the fly.
				(ind ())
				(best ()) ;; an index into the population
				(best-ind ()) ;; the best individual of the population
				(best-list ()) ;; a list of the best fitnesses
				(cycles 0)
				(ideal (funcall setup)))
;;; Create initial population
		(dotimes (p pop-size)
			(setf population (append population (funcall creator))))
			
		(dprint population "Original Population: ")
;;; Genetic Algorithm Main Loop (Algorithm 20)		
		(loop do
			(dprint "got to da top")
			;;(funcall setup)
			(setf fitnesses ())
			(dotimes (ind (length population))
				
				(dprint ind "index")
				(dprint (nth ind population) "here's the individual to evaluate")
				(setf fit (funcall evaluator (nth ind population)))
				(setf fitnesses (append fitnesses (list fit)))
				(setf *debug* t)
				(dprint fit "fitness:")
				(setf *debug* nil)
				(if (or (eql best nil) (> fit (funcall evaluator (nth best population))))
					(setf best ind)))
			(setf best-ind (nth best population))
			(setf fit (funcall evaluator best-ind))
			(dprint fitnesses "Fitnesses: ")
			(dprint (nth best fitnesses) "best:")
			(setf best-list (append best-list (list (nth best fitnesses))))
			(funcall printer population fitnesses)
			(setf q ())

;;; DAVID - do you think we can get rid of this "let" statement?
			(let ((ind1 ())
			      (ind2 ()))
				(dotimes (x  (ceiling (length population) 2))
					(setf ind1 (first (funcall selector 1 population fitnesses)))
					(setf ind2 (first (funcall selector 1 population fitnesses)))
					(setf q (append q (funcall modifier ind1 ind2)))))
			(setf population q)
			(dprint population "New Population: ")
			(setf cycles (1+ cycles))
			(print "got to da bottom")
			(print cycles)
			(print fit)
;;; DAVID - not sure how I feel about the last part of the following line ... I wonder
;;; 				if there is a better way to check for 'fit' equating to ideal?
		while (and (< cycles generations) (/= fit ideal)))
;;; DAVID - the repeated portion of Algo 20 (The GA Algo) has two (repeated) parts: looking for
;;; 				the best fit and crossover/mutation. My concern is that this function might not be
;;; 				"generic" enough. I use BEST as an index into the population, but after we mutate,
;;; 				the index is kind of useless as it might point to an element which was just modified
;;; 				(and is no longer the "best" solution. To mitigate against that, I saved off the best
;;; 				vector as BEST-IND .... but doing so may make this function less-than-generic.
;;; 				Thoughts???
;;;
;;;		(funcall printer population fitnesses) ;; problem: prints cycle AFTER best discovered
		
		(print best-ind)
		best-list
	)
)

	;;; an example way to fire up the GA.  It should easily discover
	;;; ;;; a 100% correct solution.
	;;; #|
	;;; (evolve 50 100
	;;;   :setup #'boolean-vector-sum-setup
	;;;   :creator #'boolean-vector-creator
	;;;   :selector #'tournament-selector
	;;;   :modifier #'boolean-vector-modifier
	;;;   :evaluator #'boolean-vector-evaluator
	;;;   :printer #'simple-printer)
	;;; |#

;;;;;; BOOLEAN VECTOR GENETIC ALGORTITHM

;;; Creator, Modifier, Evaluator, and Setup functions for the
;;; boolean vectors Problem.  Assume that you are evolving a bit vector
;;; *vector-length* long.  

;;; The default test function is Max-Ones.
;;;; Max-ones is defined in section 11.2.1 of "Essentials of Metaheuristics"
;;; by yours truly.  In that section are defined several other boolean functions
;;; which you should also try and report on.  Some plausible ones might
;;; include:

;;; :max-ones
;;; :trap
;;; :leading-ones
;;; :leading-ones-blocks



(defparameter *boolean-vector-length* 100)
(defparameter *boolean-problem* :max-ones)
;; perhaps other problems might include... 

(defun max-ones (vec)
	"The total number of ones in the vector vec."
	(values  (reduce #'+ vec)))

(defun leading-ones (vec)
	"Counts the number of ones in your vector,
	startin g at the beginning, until a zero is encountered."
	(dotimes (x (length vec))
		(if (equal (svref vec x) 0)
			(return x))))

(defun boolean-vector-creator ()
  "Creates a boolean-vector *boolean-vector-length* in size, filled with
random Ts and NILs, or with random 1s and 0s, your option.
(Algorithm 21.)"
	(let ((vec (make-array *boolean-vector-length* :initial-element nil)))
		(dotimes (x *boolean-vector-length*)
			(if (< (random 1.0) 0.5)
				(setf (svref vec x) 0)
				(setf (svref vec x) 1)
			)
		)
		(list vec)
	)
)

(defparameter *boolean-crossover-probability* 0.2)
(defparameter *boolean-mutation-probability* 0.01)

(defun boolean-vector-modifier (ind1 ind2)
  "Copies and modifies ind1 and ind2 by crossing them over with a uniform crossover,
then mutates the children.  *crossover-probability* is the probability that any
given allele will crossover.  *mutation-probability* is the probability that any
given allele in a child will mutate.  Mutation simply flips the bit of the allele.
(Algorithms 22 & 25.)"
	;;(dprint ind1 "ind1:")
	;;(dprint ind2 "ind2:")
	(let ((off1 (copy-seq ind1))
				(off2 (copy-seq ind2)))
	(dotimes (x (length ind1))
		;;(dprint x "Entering modifier step: ")
		(if (< (random 1.0) *boolean-crossover-probability*)
			(swap (svref off1 x) (svref off2 x))
		)
;;; These don't feel very "lispy"		
		(if (< (random 1.0) *boolean-mutation-probability*)
			(if (eql (svref off1 x) 1)
				(setf (svref off1 x) 0)
				(setf (svref off1 x) 1)))
		(if (< (random 1.0) *boolean-mutation-probability*)
			(if (eql (svref off2 x) 1)
		  	(setf (svref off2 x) 0)
		  	(setf (svref off2 x) 1)))
	)
	;;(dprint "End of modifier")
	;;(dprint ind1 "ind1:")
	;;(dprint ind2 "ind2:")
	
	(list off1 off2))
)

(defun boolean-vector-evaluator (ind1)
  "Evaluates an individual, which must be a boolean-vector, and returns
its fitness."
	(max-ones ind1)
)

(defun boolean-vector-sum-setup ()
  "Does nothing.  Perhaps you might use this to set up various
(ahem) global variables to define the problem being evaluated, I dunno."
  (let ((boolean-ideal 100))
		boolean-ideal
	)
)

;;; an example way to fire up the GA.  It should easily discover
;;; a 100% correct solution.
;;;#|
;;;(evolve 50 100
;;; 	:setup #'boolean-vector-sum-setup
;;;	:creator #'boolean-vector-creator
;;;	:selector #'tournament-selector
;;;	:modifier #'boolean-vector-modifier
;;;  :evaluator #'boolean-vector-evaluator
;;;	:printer #'simple-printer)
;;;|#




;;;;;; FLOATING-POINT VECTOR GENETIC ALGORTITHM

;;; Creator, Modifier, Evaluator, and Setup functions for the
;;; GA Max-ONES Problem.  Assume that you are evolving a vector
;;; of floating-point numbers *float-vector-length* long.  


;;; The default test function is Rastrigin.
;;;; Rastrigin is defined in section 11.2.2 of "Essentials of Metaheuristics"
;;; by yours truly.  In that section are defined several other floating-point functions
;;; which you should also try and report on.  Some plausible ones might
;;; include:

;;; :rastrigin
;;; :rosenbrock
;;; :griewank
;;; :schwefel

;;; If you were really into this, you might also try testing on
;;; rotated problems -- some of the problems above are linearly
;;; separable and rotation makes them harder and more interesting.
;;; See the end of Section 11.2.2.

;;; I have defined the minimum and maximum gene values for rastrigin
;;; as [-5.12, 5.12].  Other problems have other traditional min/max
;; ; values, consult Section 11.2.2.


(defparameter *float-vector-length* 100)

;;;(defparameter *float-problem* :rastrigin)
;;;(defparameter *float-problem* rastrigin)
(defparameter *float-min* -5.12)  ;; these will change based on the problem
(defparameter *float-max* 5.12)   ;; likewise

(defun rastrigin (a vec)
	"The <Rastrigin> function s a non-convex function used as a performance test 
problem for optimization algorithms. It is a typical example of non-linear multimodal 
function. The formula has a global minimum at x=0. Negated to account for minimum-finding
versus maximum-finding."
	(let ((n (length vec))
				(sub-total 0.0))
;;; Negated Rastrigin formula to account for *minimization* versus *maximization* problem
		(* -1.0 (+ (* a n) (reduce '+ (map 'vector #'(lambda (x) (- (* x x) (* a (cos (* 2 pi x))))) vec))))
;;;		(dotimes (x n) 
;;;			(setf sub-total (+ sub-total (- (* (svref vec x) (svref vec x)) (* a (cos (* 2 pi (svref vec x)))))))
;;;		)
;;;		(setf sub-total (* -1 (+ (* 10 n) sub-total)))
	)
)

(defun simple-sum (vec)
	"An attempt at a simple sum of the supplied vector argument." 
	(let ((n (length vec)))
		(* 1.0 (reduce '+ (map 'vector #'(lambda (x) x) vec)))
	)
)

(defun float-vector-creator ()
  "Creates a floating-point-vector *float-vector-length* in size, filled with
random numbers in the range appropriate to the given problem (Algorithm 7)."
	(let ((vec (make-array *float-vector-length* :initial-element 0.0)))
		(dotimes (x *float-vector-length*)
			(setf (svref vec x) (+ (random (- *float-max* *float-min*)) *float-min*))
		)
	(dprint (list vec) "vec:")
	)
)

(defparameter *float-crossover-probability* 0.2)
(defparameter *float-mutation-probability* 0.1)   ;; I just made up this number
(defparameter *float-mutation-variance* 0.01)     ;; I just made up this number

(defun float-vector-modifier (ind1 ind2)
  "Copies and modifies ind1 and ind2 by crossing them over with a uniform crossover,
then mutates the children.  *crossover-probability* is the probability that any
given allele will crossover.  *mutation-probability* is the probability that any
given allele in a child will mutate.  Mutation does gaussian convolution on the allele."
;;; See "Gaussian Convolution" (Algorithm 11) in the book for mutation
	;;(dprint ind1 "ind1:")
	;;(dprint ind2 "ind2:")
	(let ((off1 (copy-seq ind1))
				(off2 (copy-seq ind2))
				(n 0.0))
		(dotimes (x (length off1))
;;(dprint x "Entering modifier step: ")
			(if (< (random 1.0) *float-crossover-probability*)
				(swap (svref off1 x) (svref off2 x))
			)
;;; Alogo 11: Gaussian Convolution
			(if (< (random 1.0) *float-mutation-probability*)
				(progn 
					(loop do 
						(setf n (first (bmm 0 *float-mutation-variance*)))
;;;					(dprint n "Normal: ")
;;;					(dprint (svref off1 x) "Element1: ")
					until (and (<= *float-min* (+ n (svref off1 x))) (<= (+ n (svref off1 x)) *float-max*)))
					(setf (svref off1 x) (+ n (svref off1 x))))
			)
			(if (< (random 1.0) *float-mutation-probability*)
				(progn 
					(loop do 
						(setf n (second (bmm 0 *float-mutation-variance*)))
;;;					(dprint n "Normal: ")
;;;					(dprint (svref off2 x) "Element2: ")
					until (and (<= *float-min* (+ n (svref off2 x))) (<= (+ n (svref off2 x)) *float-max*)))
					(setf (svref off2 x) (+ n (svref off2 x))))
			)
		)
	;;(dprint "End of modifier")
	;;(dprint ind1 "ind1:")
	;;(dprint ind2 "ind2:")
	
		(list off1 off2)
	)
)

;;; Gaussian Random Sampler
(defun normal (&optional (mu 0) (sigma 1))
  (let ((u (random 1.0))
  	(v (random 1.0)))
  	    (+ mu (* sigma (sqrt (* -2 (log (- 1.0 u)))) (cos (* 2 pi v))))))

;;; Box-Muller-Marsaglia Method
(defun bmm (&optional (mu 0) (var 1))
	"Creates a pair of random numbers using the Box-Muller-Marsaglia
Method (Algorithm 12)."
	(let ((sigma (sqrt var)))
		(loop do
			(setf x (- (random 2.0) 1.0))
			(setf y (- (random 2.0) 1.0))
			(setf w (+ (* x x) (* y y)))
		until (and (> w 0) (< w 1)))
		(setf g (+ mu (* x sigma (sqrt (* -2 (float (/ (log w) w)))))))
		(setf h (+ mu (* y sigma (sqrt (* -2 (float (/ (log w) w)))))))
		(list g h)
	)
)

(defun float-vector-sum-evaluator (ind1)
  "Evaluates an individual, which must be a floating point vector, and returns
its fitness."
	(rastrigin 10.0 ind1)
;;	(simple-sum ind1)
)

(defun float-vector-sum-setup ()
  "Does nothing.  Perhaps you might use this function to set
(ahem) various global variables which define the problem being evaluated
and the floating-point ranges involved, etc.  I dunno."
  (let ((float-ideal 0))
		float-ideal
	)
)



;;; an example way to fire up the GA.  It should easily discover
;;; a 100% correct solution.
#|
(evolve 50 100
 	:setup #'float-vector-sum-setup
	:creator #'float-vector-creator
	:selector #'tournament-selector
	:modifier #'float-vector-modifier
  :evaluator #'float-vector-sum-evaluator
	:printer #'simple-printer)
|#


;;;; GP TREE CREATION CODE

;;; GP's tree individuals are considerably more complex to modify than
;;; simple vectors.  Get the GA system working right before tackling
;;; this part, trust me.



;; set up in the gp setup function -- for example, see
;; the code for gp-symbolic-regression-setup
(defparameter *nonterminal-set* nil)
(defparameter *terminal-set* nil)



;;; important hint: to use SETF to change the position of an element in a list
;;; or a tree, you need to pass into SETF an expression starting with the
;;; parent.  For example (setf (car parent) val) .... thus you HAVE to have
;;; access to the parent, not just the child.  This means that to store a
;;; "position", you have to store away the parent and the arg position of
;;; child, not the child itself.

(defun make-queue ()
  "Makes a random-queue"
  (make-array '(0) :adjustable t :fill-pointer t))
(defun enqueue (elt queue)
  "Enqueues an element in the random-queue"
  (progn (vector-push-extend elt queue) queue))
(defun queue-empty-p (queue)
  "Returns t if random-queue is empty"
  (= (length queue) 0))
(defun random-dequeue (queue)
  "Picks a random element in queue and removes and returns it.
Error generated if the queue is empty."
  (let ((index (random (length queue))))
    (swap (elt queue index) (elt queue (1- (length queue))))
    (vector-pop queue)))

(defun get-random-element (some-list)
	(dprint "random element of list ")
	(dprint some-list)
	(if (not some-list) nil
		(nth (random (length some-list)) some-list))	
)
(defun get-random-nil-element (some-list)
	(let ((value '(t)))
		(loop while (not (equal (first value)  nil)) do
			(setf value (get-random-element some-list))
		)
	value
	)
)

(defun ptc2(size)
	(if (<= size 1) (list (get-random-element *terminal-set*)) 
	(let ((temp-cell (copy-list '())) (list-of-cells (copy-list '((2)))) (current-nonterminal (copy-list '())) (new-thing (copy-list '(()))) (begining-of-list (copy-list '())))
		;; new-thing is terribly named as i wasnt in the best mood, will change before turnin
		;;it's the current cell we're modifying out of the list of cells. at first there's just a single empty cell 
		;;a "cell" is a place to put an argument
		(setf begining-of-list new-thing)
		(dotimes (n (- size 1))

			(setf current-nonterminal (get-random-element *nonterminal-set*))
			(dprint "wat, nonterminal is" current-nonterminal)
			(setf  (first new-thing) (list (first current-nonterminal)))
			(dotimes (x (first (last current-nonterminal)))
				(dprint "x is " x)
				(setf temp-cell (copy-list '()))
					
				;;the hardest step, figureing out how to add this thing to the end of the darn linked list without screwing things up
				(push (copy-list '()) temp-cell)
				(setf (cdr (last list-of-cells)) (list temp-cell))	
				(setf (cdr (last (first new-thing))) temp-cell)
				
			)
		
			(dprint "list of cells" list-of-cells)
			(dprint "first of newthing" (first new-thing))
		
			(setf new-thing (get-random-nil-element list-of-cells))
			;;(setf (first new-thing) (random 400))
			
			(dprint "HEY new-thing"  new-thing)
			;;(setf new-thing (first list-of-cells))
			(dprint begining-of-list "begining of list")
		)
		(loop for x in list-of-cells do
			(if (equal (first x) nil) (setf (first x) (list (get-random-element *terminal-set*)))))
		(first begining-of-list)
	))
)



(defparameter *size-limit* 20)
(defun gp-creator ()
  "Picks a random number from 1 to 20, then uses ptc2 to create
a tree of that size"
	(print "HEY HERE'S SPMETJNG I MADE")
	(print (list (ptc2 (+ (random *size-limit*) 1))))

)



;;; GP TREE MODIFICATION CODE

(defun num-nodes (tree)
  "Returns the number of nodes in tree, including the root"
    (let ((sum 0));;nice simple recursive implementation for counting of nodes
	(dotimes (n (length tree))
    		(if (atom (nth n tree)) (incf sum)  (setf sum (+ sum (num-nodes (nth n tree)))))
	)
	sum;;return value	
    )
)

(defun pre-process (tree)
    (setf tree (macroexpand tree))
    (let ((sum 0));;nice simple recursive implementation for counting of nodes
	(dotimes (n (length tree))
    		(if (atom (nth n tree)) (incf sum)  (progn (setf (nth n tree) (macroexpand (nth n tree))) (setf (nth n tree) (pre-process (nth n tree)))))
	)
	tree;;return value	
    )
)

(defun nth-subtree-parent (tree n)
  "Given a tree, finds the nth node by depth-first search though
the tree, not including the root node of the tree (0-indexed). If the
nth node is NODE, let the parent node of NODE is PARENT,
and NODE is the ith child of PARENT (starting with 0),
then return a list of the form (PARENT i).  For example, in
the tree (a (b c d) (f (g (h) i) j)), let's say that (g (h) i)
is the chosen node.  Then we return ((f (g (h) i) j) 0).

If n is bigger than the number of nodes in the tree
 (not including the root), then we return n - nodes_in_tree
 (except for root)."

  ;;; this is best described with an example:
  ;    (dotimes (x 12)
  ;           (print (nth-subtree-parent
  ;                       '(a (b c) (d e (f (g h i j)) k))
  ;                        x)))
  ;;; result:
  ;((A (B C) (D E (F (G H I J)) K)) 0) 
  ;((B C) 0) 
  ;((A (B C) (D E (F (G H I J)) K)) 1) 
  ;((D E (F (G H I J)) K) 0) 
  ;((D E (F (G H I J)) K) 1) 
  ;((F (G H I J)) 0) 
  ;((G H I J) 0) 
  ;((G H I J) 1) 
  ;((G H I J) 2) 
  ;((D E (F (G H I J)) K) 2) 
  ;0 
  ;1 
  ;NIL
;;	(print "recursive call on")
;;	(print n)
	
	
	;;NOTE: not everything here has a purpose, and there are some really stupid looking !@#$ going on. 
	;;this is a result of me changing the algorithm to solve this as i was debugging it
	;;There are left over variables for trying to make ideas for how to actually 
	;;solve this problem i'll clean it up when we have overhead to do so.
	;;inspired off of num-nodes, and uses num nodes. 
	;; will only do a recursive call iff that subtree is big enough that it MUST have the desired subtree in it
    (if (>= (+ n 1) (num-nodes tree)) (return-from nth-subtree-parent (+ (- n (num-nodes tree)) 1 ))) 
    (let ((sum 0) (k 0) (temp 0) (found '()) (parent '()));;nice simple recursive implementation for counting of nodes
        (dotimes (my-n (length tree))
                (setf k my-n)
		(if (= n -1) (return-from nth-subtree-parent (list tree (- my-n 1) )))
		(if (atom (nth my-n tree)) (progn (decf n) )   
			(progn 
				(dprint "before second return")
				(if (= n -1) (return-from nth-subtree-parent (list (nth my-n tree) 0)))
				(dprint "after second return")
				(setf temp n)
				(setf sum (+ sum (num-nodes (nth my-n tree))))
				
				(dprint sum)
				(dprint n)
				(dprint "before third return")
				(if (> sum (+ n 1)) (return-from nth-subtree-parent (nth-subtree-parent (nth my-n tree) (+ n 0))) ()) 
				(dprint "after third return")
				(setf n (- n sum))
				(setf sum 0)
			))
		        
	)
     )
    ;;; IMPLEMENT ME

  )


(defparameter *mutation-size-limit* 10)
(defun crossover-gp (ind1 ind2)
	(dprint "before")
	
	(eval ind1)
	(eval ind2)
	(dprint "ind1: ")
	(dprint ind1)
	(dprint "ind2:")
	(dprint ind2)
	(dprint "after first eval")
	(let (  (first-index 0)
		(second-index 0)
		(new-tree1 '())
		(new-tree2 '())
		(subtree1 '())
		(subtree2 '()))
	
		(if (> (length ind1) 1)	
			(setf subtree1 (nth-subtree-parent ind1 (random (- (num-nodes ind1) 1))))
			(setf subtree1 (list (list ind1) -1)))
		(if (> (length ind2) 1)
			(setf subtree2 (nth-subtree-parent ind2 (random (- (num-nodes ind2) 1))))
			(setf subtree2 (list (list ind2) -1)))
		
		(setf first-index (+ (second subtree1) 1))
		(setf second-index (+ (second subtree2) 1))

		(setf new-tree1 (copy-tree (nth first-index (first subtree1))))
		(setf new-tree2 (copy-tree (nth second-index (first subtree2))))
		(setf (nth first-index (first subtree1)) new-tree2) 
		(setf (nth second-index (first subtree2)) new-tree1) 
		(if (and (= 1 (length ind1) ) (not (atom (first ind1)))) (setf ind1 (first ind1)))
		(if (and (= 1 (length ind2) ) (not (atom (first ind1)))) (setf ind2 (first ind2)))
		(list ind1 ind2)
	) 
	(dprint "beofre2")
	(dprint "ind1: ")
        (dprint ind1)
        (dprint "ind2:")
        (dprint ind2)
	
	(eval ind1)
	(eval ind2)
	
	(dprint "after2")
	(list ind1 ind2)
)
(defun modify-tree (ind1)
	(dprint "BEFORE")
	(dprint ind1)
	(dprint "hello world")
	(dprint (num-nodes ind1))
	(let (  (first-index 0)
		(new-tree (ptc2 (+ (random *mutation-size-limit*) 1)))
		(subtree1 '()))
		(if (> (length ind1) 1) 
			(setf subtree1 (nth-subtree-parent ind1 (random (- (num-nodes ind1) 1))))
			(return-from modify-tree new-tree))
		;;(print "hi1")
		;;(print "nth-subtree returned")
		;;(print subtree1)
		(setf first-index (+ (second subtree1) 1))
		;;(print "hi2")
		(setf (nth first-index (first subtree1)) new-tree) 
		;;(print "hi3")
		(dprint "AFTER")
		(dprint ind1)
		(if (atom ind1) (print (/ 1 0)))
		ind1			
	)
	(dprint "before eval")
	(eval ind1)
	(dprint "after eval")
	ind1
)
(defun gp-modifier (ind1 ind2)
  "Flips a coin.  If it's heads, then ind1 and ind2 are
crossed over using subtree crossover.  If it's tails, then
ind1 and ind2 are each mutated using subtree mutation, where
the size of the newly-generated subtrees is pickedc at random
from 1 to 10 inclusive.  Doesn't damage ind1 or ind2.  Returns
the two modified versions as a list."
	(setf ind1 (copy-tree ind1))
	(setf ind2 (copy-tree ind2))
	(if (= (random 2) 0)
		(crossover-gp ind1 ind2)
		(list (modify-tree ind1) (modify-tree ind2))
	
	)	
    ;;; IMPLEMENT ME
)






;;; SYMBOLIC REGRESSION
;;; This problem domain is similar, more or less, to the GP example in
;;; the lecture notes.  Your goal is to make a symbolic expression which
;;; represents a mathematical function consisting of SIN COS, EXP,
;;; +, -, *, and % (a version of / which doesn't give an error on
;;; divide-by-zero).  And also the function X, which returns a current
;;; value for the X variable.
;;;
;;; In symbolic regression, we generate 20 (x,y) pairs produced at
;;; random which fit the expression y = x^4 + x^3 + x^2 + x.  You can
;;; make up another function is you like, but that's the one we're going
;;; with here.  Using just these data pairs, the GP system will attempt
;;; to ascertain the function.  It will generate possible expressions
;;; as its individuals, and the fitness of an expression is how closely,
;;; for each X value, it gives the correct corresponding Y value.
;;;
;;; This is called SYMBOLIC regression because we're actually learning
;;; the mathematical expression itself, including transcendental functions
;;; like SIN and COS and E^.  As opposed to statistical linear or
;;; quadratic curve-fitting regressions which just try to learn the
;;; linear or quadratic parameters etc.
;;;
;;; An example 100% optimal solution:
;;;
;;; (+ (* (x) (* (+ (x) (* (x) (x))) (x))) (* (+ (x) (cos (- (x) (x)))) (x)))



;;; GP SYMBOLIC REGRESSION SETUP
;;; (I provide this for you)
(defun nothing (one two three)
	one
)
(defparameter *num-vals* 20)
(defparameter *vals* nil) ;; gets set in gp-setup

(defun gp-symbolic-regression-setup ()
  "Defines the function sets, and sets up vals"

  (setq *nonterminal-set* '((+ 2) (- 2) (* 2) (% 2) (sin 1) (cos 1) (exp 1) (nothing 3) ))
  (setq *terminal-set* '(x))

  (setq *vals* nil)
  (dotimes (v *num-vals*)
    (push (1- (random 2.0)) *vals*))
	1
)

(defun poly-to-learn (x) (+ (* x x x x) (* x x x) (* x x) x))

;; define the function set
(defparameter *x* nil) ;; to be set in gp-evaluator
(defun x () *x*)
(defun % (x y) (if (= y 0) 0 (/ x y)))  ;; "protected division"
;;; the rest of the functions are standard Lisp functions




;;; GP SYMBOLIC REGRESSION EVALUATION

(defun gp-symbolic-regression-evaluator (ind)
  "Evaluates an individual by setting *x* to each of the
elements in *vals* in turn, then running the individual and
get the output minus (poly-to-learn *x*).  Take the
absolute value of the this difference.  The sum of all such
absolute values over all *vals* is the 'raw-fitness' Z.  From
this we compute the individual's fitness as 1 / (1 + z) -- thus
large values of Z are low fitness.  Return the final
individual's fitness.  During evaluation, the expressions
evaluated may overflow or underflow, or produce NaN.  Handle all
such math errors by
returning most-positive-fixnum as the output of that expression."
  ;;; hint:
  ;;; (handler-case
  ;;;  ....
  ;;;  (error (condition)
  ;;;     (format t "~%Warning, ~a" condition) most-positive-fixnum))
	(handler-case

	(let ((sum 0.0))
		(loop for x in *vals* do (progn
			(dprint "some value in vals")
			(dprint x)
			(dprint "ind is ")
			(dprint ind)
			(setf *x* x)
				(setf sum 
					(+ sum (abs (- (float (eval ind)) (poly-to-learn x))))))	
		)
		(if (not (equal sum nil)) (/ 1.0 (+ sum 1)) (0))
	)
	(error (condition)
         		   (format t "~%Warning, ~a" condition) (return-from gp-symbolic-regression-evaluator 0)))
 
	 ;;; IMPLEMENT ME

  )


;;; Example run
#|
(evolve 50 500
 	:setup #'gp-symbolic-regression-setup
	:creator #'gp-creator
	:selector #'tournament-selector
	:modifier #'gp-modifier
  :evaluator #'gp-symbolic-regression-evaluator
	:printer #'simple-printer)
|#





;;; GP ARTIFICIAL ANT CODE
;;; for this part you'll need to implement both the evaluator AND
;;; make up the function that form the function set.

;;; In the artificial ant problem domain, you'll be evolving an s-expression
;;; which moves an ant around a toroidal map shown below.  The ant starts out
;;; at (0,0), facing east (to the right).  The functions in
;;; the expression are:
;;; (if-food-ahead --then-- --else--)   If food is directly ahead of the ant,
;;;                                     then evaluate the THEN expression, else
;;;                                     evaluate the ELSE expression
;;; (progn2 --item1-- --item2--)        Evaluate item1, then evaluate item 2
;;; (progn3 --item1-- --item2-- --item3--)  Evaluate item1, then item2, then item3
;;; (move)                              Move forward one unit
;;;                                     If you pass over a food pellet, it is eaten
;;;                                     and removed from the map.
;;; (left)                              Turn left
;;; (right)                             Turn right
;;;
;;;
;;; the way a tree is evaluated is as follows.  The whole tree is evaluated,
;;; and the functions move the ant about.  If the ant has not made a total of
;;; 600 MOVE, LEFT, and RIGHT operations yet, then the tree is evaluated AGAIN
;;; moving the ant some more from its current position.  This goes on until
;;; 600 operations have been completed, at which time the ant will not move
;;; any more.  If 600 operations are completed in the middle of the evaluation
;;; of a tree, the simplest approach is to have the MOVE command simply
;;; "stop working" so the ant doesn't gobble up any more pellets accidentally.

;;; The fitness of the artificial ant is the number of pellets it ate in the
;;; 600-operation period.  Maps are reset between evaluations of different
;;; individuals of course.

;;; Here's an optimal ant program (there are many such):
;;  (progn3 (if-food-ahead (move) (progn2 (left) (progn2 (left)
;;             (if-food-ahead (move) (right))))) (move) (right)))

;;; This is a hard problem for GP and you may need to run many times before you
;;; get a 100% perfect answer.

;;; Note that in my thesis it says 400 moves and not 600.  We're going with
;;; 600 here.  It's easier.



;;; our ant's food trail map
(defparameter *map-strs* '(
".###............................"
"...#............................"
"...#.....................###...."
"...#....................#....#.."
"...#....................#....#.."
"...####.#####........##........."
"............#................#.."
"............#.......#..........."
"............#.......#..........."
"............#.......#........#.."
"....................#..........."
"............#..................."
"............#................#.."
"............#.......#..........."
"............#.......#.....###..."
".................#.....#........"
"................................"
"............#..................."
"............#...#.......#......."
"............#...#..........#...."
"............#...#..............."
"............#...#..............."
"............#.............#....."
"............#..........#........"
"...##..#####....#..............."
".#..............#..............."
".#..............#..............."
".#......#######................."
".#.....#........................"
".......#........................"
"..####.........................."
"................................"))

(defparameter *map-height* 32)
(defparameter *map-width* 32)


;;; some useful functions for you

(defun make-map (lis)
  "Makes a map out of a string-list such as *MAP-STRS*"
  (let ((map (make-array (list (length (first lis)) (length lis)))))
    (dotimes (y (length lis) map)
      (dotimes (x (length (elt lis y)))
	(setf (aref map x y)
	      (cond ((equalp #\# (elt (elt lis y) x)) nil)
		    (t t)))))))

(defun direction-to-arrow (dir)
  "Returns a character which represents a given direction -- might
be useful for showing the movement along a path perhaps..."
  (cond ((= dir *n*) #\^)
	((= dir *s*) #\v)
	((= dir *e*) #\>)
	(t #\<)))

(defun maparray (function array &optional (same-type nil))
  "Maps function over array in its major order.  If SAME-TYPE,
then the new array will have the same element type as the old
array; but if a function returns an invalid element then an error
may occur.  If SAME-TYPE is NIL, then the array will accommodate
any type."
  (let ((new-array (apply #'make-array (array-dimensions array)
			  :element-type (if same-type
					    (array-element-type array)
					  t)
			  :adjustable (adjustable-array-p array)
			  (if (vectorp array)
			      `(:fill-pointer ,(fill-pointer array))
			    nil))))
    (dotimes (x (array-total-size array) new-array)
      (setf (row-major-aref new-array x)
	    (funcall function (row-major-aref array x))))))

(defun print-map (map)
  "Prints a map, which must be a 2D array.  If a value in the map
is T (indicating a space), then a '.' is printed.  If a value in the map
is NIL (indicating a food pellet), then a '#' is printed.  If a value in
the map is anything else, then it is simply PRINTed.  This allows you to
consider spaces to be all sorts of stuff in case you'd like to print a
trail of spaces on the map for example.  Returns NIL."
  (let ((dim (array-dimensions map)))
    (dotimes (y (second dim) nil)
      (format t "~%")
      (dotimes (x (first dim))
	(format t "~a"
		(let ((v (aref map x y)))
		  (cond ((equal v t) #\.)
			((null v) #\#)
			(t v))))))))


;; The four directions.  For relative direction, you might
;; assume that the ant always PERCEIVES things as if it were
;; facing north.
(defconstant *n* 0)
(defconstant *e* 1)
(defconstant *s* 2)
(defconstant *w* 3)

(defmacro absolute-direction (relative-dir ant-dir)
  "If the ant is facing ANT-DIR, then converts the perceived
RELATIVE-DIR direction into an absolute ('true') direction
and returns that."
  `(mod (+ ,relative-dir ,ant-dir) 4))

(defmacro x-pos-at (x-pos absolute-dir &optional (steps 1))
  "Returns the new x position if one moved STEPS steps the absolute-dir
direction from the given x position.  Toroidal."

   `(mod (cond ((= (mod ,absolute-dir 2) *n*) ,x-pos)         ;; n or s
	      ((= ,absolute-dir *e*) (+ ,x-pos ,steps))     ;; e
	      (t (+ ,x-pos (- ,steps) *map-width*)))         ;; w
	*map-width*))

(defmacro y-pos-at (y-pos absolute-dir &optional (steps 1))
  "Returns the new y position if onee moved STEPS steps in the absolute-dir
direction from the given y position.  Toroidal."
  `(mod (cond ((= (mod ,absolute-dir 2) *e*) ,y-pos)        ;; e or w
	      ((= ,absolute-dir *s*) (+ ,y-pos ,steps))     ;; s
	      (t (+ ,y-pos (- ,steps) *map-height*)))       ;; n
	*map-height*))


(defparameter *current-move* 0 "The move # that the ant is at right now")
(defparameter *num-moves* 600 "How many moves the ant may make")
(defparameter *current-x-pos* 0 "The current X position of the ant")
(defparameter *current-y-pos* 0 "The current Y position of the ant")
(defparameter *current-ant-dir* *e* "The current direction the ant is facing")
(defparameter *eaten-pellets* 0 "How many pellets the ant has eaten so far")
(defparameter *map* (make-map *map-strs*) "The ant's map")




;;; the function set you have to implement

(defmacro if-food-ahead (then else)
  "If there is food directly ahead of the ant, then THEN is evaluated,else ELSE is evaluated"
  ;; because this is an if/then statement, it MUST be implemented as a macro.
   `(if 
      (aref *map* (x-pos-at *current-x-pos* *current-ant-dir*) (y-pos-at *current-y-pos* *current-ant-dir*))
      ,else
	  ,then)

    
    ;;; IMPLEMENT ME
)

(defun progn2 (arg1 arg2)
    "Evaluates arg1 and arg2 in succession, then returns the value of arg2"
    (declaim (ignore arg1))
    arg2)  ;; ...though in artificial ant, the return value isn't used ... 

(defun progn3 (arg1 arg2 arg3)
  "Evaluates arg1, arg2, and arg3 in succession, then returns the value of arg3"
  (declaim (ignore arg1 arg2))
  arg3)  ;; ...though in artificial ant, the return value isn't used ...

(defun move ()
  "If the move count does not exceed *num-moves*, increments the move count
and moves the ant forward, consuming any pellet under the new square where the
ant is now.  Perhaps it might be nice to leave a little trail in the map showing
where the ant had gone."
		(incf *current-move*)
		(if (< *current-move* *num-moves*) (progn
			
			
			;;mark my path
			(setf (aref *map* *current-x-pos* *current-y-pos*) (direction-to-arrow *current-ant-dir* ))
			;;actually move
			(setf *current-x-pos* (x-pos-at *current-x-pos* *current-ant-dir*))
			(setf *current-y-pos* (y-pos-at *current-y-pos* *current-ant-dir*))
			;;did i find food?
			(if (equal (aref *map* *current-x-pos* *current-y-pos*) nil)
				(incf *eaten-pellets*)
				
			)
		))
      ;;; IMPLEMENT ME
  )


(defun left ()
  "Increments the move count, and turns the ant left"
	(incf *current-move*)
	(if (< *current-move* *num-moves*)
		(setf *current-ant-dir* (absolute-direction -1 *current-ant-dir*))
	)
		
      ;;; IMPLEMENT ME
)

(defun right ()
  "Increments the move count, and turns the ant right"
	(incf *current-move*)
	(if (< *current-move* *num-moves*)
		(setf *current-ant-dir* (absolute-direction 1 *current-ant-dir*))
	)
      ;;; IMPLEMENT ME
)

(defparameter *nonterminal-set* nil)
(defparameter *terminal-set* nil)

;; I provide this for you
(defun gp-artificial-ant-setup ()
  "Sets up vals"
  (setq *nonterminal-set* '((if-food-ahead 2) (progn2 2) (progn3 3) ))
  (setq *terminal-set* '(left right move))
  (setq *map* (make-map *map-strs*))
  (setq *current-move* 0)
  (setq *eaten-pellets* 0)
  (setf *num-moves 600)
  (setf *size-limit* 20)
  89
)


;; you'll need to implement this as well

(defun gp-artificial-ant-evaluator (ind)
  "Evaluates an individual by putting it in a fresh map and letting it runfor *num-moves* moves.  The fitness is the number of pellets eaten -- thusmore pellets, higher (better) fitness."
 	(setf *map* (make-map *map-strs*))
  	(setf *current-move* 0)
  	(setf *eaten-pellets* 0)	
	(setf *current-ant-dir* 0)
	(setf *current-x-pos* 0)
	(setf *current-y-pos* 0)
	(dprint "before")
	(dprint ind)
	(setf ind (copy-tree ind))
	(dotimes (x 10)
	(setf ind (pre-process ind)))
	(dprint "preprocess 1")
	(dprint ind)
	(setf ind (pre-process ind))
	(dprint "preprocess 2")
	(dprint ind)
	(block dotimes
		(dotimes (x *num-moves*)
			(eval ind)
			(if (> *current-move* *num-moves*) (return-from dotimes nil))
	))
	*eaten-pellets*
	
      ;;; IMPLEMENT ME
)

;; you might choose to write your own printer, which prints out the best
;; individual's map.  But it's not required.

#|
(evolve 50 500
 	:setup #'gp-artificial-ant-setup
	:creator #'gp-creator
	:selector #'tournament-selector
	:modifier #'gp-modifier
  :evaluator #'gp-artificial-ant-evaluator
	:printer #'simple-printer)
|#

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Test Code for boolean-vector-creator
(defun vec-test ()
	(let ((vec-1 (boolean-vector-creator))
				(vec-2 (boolean-vector-creator))
				(new-vec ()))
		(dprint vec-1 "Vector 1")
		(dprint (boolean-vector-evaluator vec-1) "Fitness 1")
		(dprint vec-2 "vector 2")
		(dprint (boolean-vector-evaluator vec-2) "Fitness 2")
		(setf new-vec (boolean-vector-modifier vec-1 vec-2))
		new-vec
;;;		(dprint new-vec "Modified ")
	)
)

;;;(vec-test)
;;;
(defun e-test ()
(setf *debug* nil)
  (evolve 50 100
		:setup #'boolean-vector-sum-setup
 		:creator #'boolean-vector-creator
 		:selector #'tournament-selector
 		:modifier #'boolean-vector-modifier
 		:evaluator #'boolean-vector-evaluator
 		:printer #'simple-printer)
)

; Some simple examples/unit tests of the Rastrigin function.
 (assert (= 0.0 (print (rastrigin 10 '(0 0))))) ; This is the global minimum.
 (assert (= 0.0 (rastrigin 10 '(0.0 0.0)))) ; This is the global minimum.
 (assert (= 0.0 (rastrigin 10.0 '(0.0 0.0)))) ; This is the global minimum.
 (assert (= 0.0 (rastrigin 10 '(0 0 0 0 0 0)))) ; This is the global minimum.

(defun f-test ()
	(setf *debug* t)
	(evolve 50 100
 		:setup #'float-vector-sum-setup
		:creator #'float-vector-creator
		:selector #'tournament-selector
		:modifier #'float-vector-modifier
  	:evaluator #'float-vector-sum-evaluator
		:printer #'simple-printer)
)
(defun ptc2-test ()
	(setf *terminal-set* '(x))
	(setf *nonterminal-set* '((+ 2) (- 2) (* 2) (% 2) (sin 1) (cos 1) (exp 1)))
 

	(print "here's (ptc2 1)")
	(print (ptc2 1))
	(print "here's (ptc2 10)")
	(print (ptc2 10))

)
(defun ptc2-test-ant ()
	(gp-artificial-ant-setup)
 

	(print "here's (ptc2 1)")
	(print (ptc2 1))
	(print "here's (ptc2 10)")
	(print (ptc2 10))

	(dotimes (n 10000)
		(eval (ptc2 (+ 1 (random 20)))))
)
(defun num-nodes-test ()
	(setf expression (ptc2 10))
	(print "expression:")
	(print expression)
	(print "num nodes is:")
	(print (num-nodes expression))
)
(defun test-subtree ()
      (let ((my-x 0))
	
	(dotimes (my-x 12) 
             (print (nth-subtree-parent '(a (b) (1 2) (d e (f (g h i j)) k)) my-x ))))
) 
(defun crossover-test ()
	(dotimes (blah 10)
		(print "origional:")
		(print '(1 (2 (5 (6 7)))))
		(print '(a b c (d e (f g h))))
		(print (crossover-gp (copy-tree '(1 (2 (5 (6 7))))) (copy-tree '(a b c (d e (f g h)))))))
	
)

;;(defun poly-to-learn (x) (+ (* x x x x) (* x x x) (* x x) x))

;; define the function set
;;(defparameter *x* nil) ;; to be set in gp-evaluator
;;(defun x () *x*)
;;(defun % (x y) (if (= y 0) 0 (/ x y)))  ;; "protected division"
;;; the rest of the functions are standard Lisp functions




;;; GP SYMBOLIC REGRESSION EVALUATION

;;(defun gp-symbolic-regression-evaluator (ind)
 
(defun modify-tree-test ()
	(setf *debug* nil)
	(setf *terminal-set* '(x))
	(setf *nonterminal-set* '((+ 2) (- 2) (* 2) (% 2) (sin 1) (cos 1) (exp 1)))
 	(dprint "doing that mutation 100 times")
	(dotimes (blah 100)
		(print (modify-tree (copy-tree '(x)))))	
	(setf *random-tree* (ptc2 5))
	(dprint "randomtere:")
	(dprint *random-tree*)
	(dotimes (blah 10)
		(dprint "origional:")
		(dprint *random-tree*)
		(dprint "modified:")
		(dprint (modify-tree *random-tree*))
))
(defun test-gp-evaluator ()
	(gp-symbolic-regression-setup)	
	(print "error with function X")
	(print (gp-symbolic-regression-evaluator '(X)))	
	(print "error with perfect function")
	(print (gp-symbolic-regression-evaluator '(+ (* (x) (x) (x) (x)) (* (x) (x) (x)) (* (x) (x)) (x))))

)
(defun test-gp-symbolic ()
(evolve 50 500
 	:setup #'gp-symbolic-regression-setup
	:creator #'gp-creator
	:selector #'tournament-selector
	:modifier #'gp-modifier
  :evaluator #'gp-symbolic-regression-evaluator
	:printer #'simple-printer)
	(gp-symbolic-regression-setup)	
)
(defun test-float ()
(evolve 100 10000
 	:setup #'float-vector-sum-setup
	:creator #'float-vector-creator
	:selector #'tournament-selector
	:modifier #'float-vector-modifier
  :evaluator #'float-vector-sum-evaluator
	:printer #'simple-printer))
(setf *debug* nil)

(defun test-ant ()
	(evolve 10 50
 	:setup #'gp-artificial-ant-setup
	:creator #'gp-creator
	:selector #'tournament-selector
	:modifier #'gp-modifier
  :evaluator #'gp-artificial-ant-evaluator
	:printer #'simple-printer)

)
;;     (sb-sprof:with-profiling (:max-samples 1000000
;;                               :mode :alloc
;;                               :report :flat)
;;       (test-ant))
(test-ant)
;;(test-gp-symbolic);;
;;(ptc2-test)
;;(num-nodes-test)
;;(test-subtree)
;;(crossover-test)
;;(modify-tree-test)
;;(test-gp-evaluator)
;;(test-ant)
;;(ptc2-test-ant)
(print-map *map*)
;;(print "hey i finished at least")
