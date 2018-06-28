/*----------------------------------------------------------------------------

  LSD - Line Segment Detector on digital images

  This code is part of the following publication and was subject
  to peer review:

    "LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
    Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
    Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
    http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd

  Copyright (c) 2007-2011 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** @file lsd.c
    
    LSD module code
    @author rafael grompone von gioi <grompone@gmail.com>

    Java port 
    @author chris <anfractuosity@anfractuosity.com>

 */
/*----------------------------------------------------------------------------*/

/** @mainpage LSD code documentation

    This is an implementation of the Line Segment Detector described
    in the paper:

      "LSD: A Fast Line Segment Detector with a False Detection Control"
      by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
      and Gregory Randall, IEEE Transactions on Pattern Analysis and
      Machine Intelligence, vol. 32, no. 4, pp. 722-732, April, 2010.

    and in more details in the CMLA Technical Report:

      "LSD: A Line Segment Detector, Technical Report",
      by Rafael Grompone von Gioi, Jeremie Jakubowicz, Jean-Michel Morel,
      Gregory Randall, CMLA, ENS Cachan, 2010.

    The version implemented here includes some further improvements
    described in the following publication, of which this code is part:

      "LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
      Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
      Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
      http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd

    The module's main function is lsd().

    The source code is contained in two files: lsd.h and lsd.c.

    HISTORY:
    - version 1.6 - nov 2011:
                              - changes in the interface,
                              - max_grad parameter removed,
                              - the factor 11 was added to the number of test
                                to consider the different precision values
                                tested,
                              - a minor bug corrected in the gradient sorting
                                code,
                              - the algorithm now also returns p and log_nfa
                                for each detection,
                              - a minor bug was corrected in the image scaling,
                              - the angle comparison in "isaligned" changed
                                from < to <=,
                              - "eps" variable renamed "log_eps",
                              - "lsd_scale_region" interface was added,
                              - minor changes to comments.
    - version 1.5 - dec 2010: Changes in 'refine', -W option added,
                              and more comments added.
    - version 1.4 - jul 2010: lsd_scale interface added and doxygen doc.
    - version 1.3 - feb 2010: Multiple bug correction and improved code.
    - version 1.2 - dec 2009: First full Ansi C Language version.
    - version 1.1 - sep 2009: Systematic subsampling to scale 0.8 and
                              correction to partially handle "angle problem".
    - version 1.0 - jan 2009: First complete Megawave2 and Ansi C Language
                              version.

    @author rafael grompone von gioi <grompone@gmail.com>
 */
/*----------------------------------------------------------------------------*/
package LSD;
import java.util.ArrayList;
import java.awt.Point;
class rect_itr {
	double[] vx; /* 以圆形顺序矩形的角X坐标 */
	double[] vy; /* 矩形的Y角坐标按循环顺序排列 */
	double ys, ye; /* 开始和结束当前'列'的Y值 */
	int x, y; /* 当前探索像素的坐标 */

	rect_itr() {
		vx = new double[4];
		vy = new double[4];
	}
}

class image_double {
	double[] data;
	int xsize, ysize;

	image_double(int xsize, int ysize) {//像素和分辨率
		this.xsize = xsize;
		this.ysize = ysize;

		data = new double[xsize * ysize];
	}

	image_double(int xsize, int ysize, double[] data) {
		this.xsize = xsize;
		this.ysize = ysize;

		this.data = data;
	}

}

/*----------------------------------------------------------------------------*/
/**
 * 链接的坐标列表
 */
class coorlist {
	int x, y;
	coorlist next;
};

public class LSD {
	
	double[] inv = new double[TABSIZE]; /*
	 * table to keep computed inverse
	 * values
	 */
	
	/** ln(10) */
	double M_LN10 = 2.30258509299404568402;

	
	
	/** PI */
	double M_PI = 3.14159265358979323846;

	int FALSE = 0;

	int TRUE = 1;

	/** Label for pixels with undefined gradient. */
	double NOTDEF = -1024.0;

	/** 3/2 pi */
	double M_3_2_PI = 4.71238898038;

	/** 2 pi */
	double M_2__PI = 6.28318530718;

	/** Label for pixels not used in yet. */
	int NOTUSED = 0;

	/** Label for pixels already used in detection. */
	int USED = 1;

	/*----------------------------------------------------------------------------*/
	/**
	 * Chained list of coordinates.
	 */
	class Coorlist {
		int x, y;
		Coorlist next;
	};

	void error(String msg) {
		System.err.println("LSD Error: " + msg);
		throw new RuntimeException("major error");
	}

	/*----------------------------------------------------------------------------*/
	/**
	 *加倍相对误差因素
	 */
	double RELATIVE_ERROR_FACTOR = 100.0;

	private double calculateMachineEpsilonDouble() {
		float machEps = 1.0f;

		do
			machEps /= 2.0f;
		while ((double) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Compare doubles by relative error.
  *浮点运算后产生的舍入误差取决于
  *完成的具体操作。 相同的数字通过不同的计算
  *算法可能会出现不同的舍入错误。 对于一个有用的
  *比较，相对舍入误差的估计应该是
  *考虑并将其与EPS乘以因子时间进行比较。 因素应该是
  *与计算链中的累计舍入误差有关。
  *这里，作为简化，使用固定因子。
  */
	
	boolean double_equal(double a, double b) {
		double abs_diff, aa, bb, abs_max;

		/* trivial case */
		if (a == b)
			return true;

		abs_diff = Math.abs(a - b);
		aa = Math.abs(a);
		bb = Math.abs(b);
		abs_max = aa > bb ? aa : bb;

		/*
		* DBL_MIN是最小的归一化数字，因此是最小的数字
        *其相对误差受DBL_EPSILON限制。 对于较小的数字，
        *使用与DBL_MIN相同的量化步骤。 然后，为
        *较小的数字，应该通过计算一个有意义的“相对”误差
        *用DBL_MIN除以差值。
		 */
		if (abs_max < -Double.MAX_VALUE)
			abs_max = -Double.MAX_VALUE;

		/*  如果相对误差<=因子x eps，则相等 */
		return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * calculateMachineEpsilonDouble());
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * 计算点（x1，y1）和点（x2，y2）之间的欧几里得距离。
	 */
	static double dist(double x1, double y1, double x2, double y2) {
		return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	}

	/*----------------------------------------------------------------------------*/
	/*----------------------- 'n元组'列表的数据类型 ------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	*'n元组'数据类型列表
	*
	* n元组列表“ntl”的第j个n元组的第i个元素是
	*访问：
	*
	* ntl.values [i + j * ntl.dim]
	*
	* n元组（n）的维数是：
	*
	* ntl.dim
	*
	*列表中的n元组数量为：
	*
	* ntl.size
	*
	*可用列表存储在列表中的n元组的最大数量
	*给定时间分配的内存由下式给出：
	*
	* ntl.max_size
	*/
	class ntuple_list {
		int size;
		int max_size;
		int dim;
		double[] values;

		/*----------------------------------------------------------------------------*/
		/**
		 *创建一个n元组列表并为一个元素分配内存。
		 * 
		 * @param dim
		 *           n元组的维度（n）。
		 */
		ntuple_list(int dim) {

			/* check parameters */
			if (dim == 0)
				error("new_ntuple_list: 'dim' must be positive.");

			/* initialize list */
			size = 0;
			max_size = 1;
			this.dim = dim;

		     // *获取元组的内存
			// n_tuple.values = new ArrayList（）; （double *）malloc（*
			// sizeof（double）
			//）;
			//如果（n_tuple.values == NULL）错误（“内存不足”）;

			values = new double[dim * max_size];

		}

	}

	/*----------------------------------------------------------------------------*/
	/**
	 * 放大n元组列表的分配内存。
	 */
	void enlarge_ntuple_list(ntuple_list n_tuple) {
		 

	//	/ *检查参数* /
		// if（n_tuple == null || n_tuple.values == null || n_tuple.max_size ==
		// 0）
		// error（“enlarge_ntuple_list：invalid n-tuple。”）;

		/// *元组的重复数量* /
		n_tuple.max_size *= 2;

		// *重新分配内存* /

		//System.out.println("THIS IS ACTUALLY WRONG!!!!!!!!!!");
		int oldlen = n_tuple.values.length;
		
		
		double [] arr = new double[n_tuple.dim * n_tuple.max_size];
		
		
		for(int i=0; i < oldlen; i++){
			arr[i] = n_tuple.values[i];
		}

		n_tuple.values = arr;
		
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * 将一个7元组添加到n元组列表中。
	 */
	void add_7tuple(ntuple_list out, double v1, double v2, double v3,
			double v4, double v5, double v6, double v7) {
		/* check parameters */
		if (out == null)
			error("add_7tuple: invalid n-tuple input.");
		if (out.dim != 7)
			error("add_7tuple: the n-tuple must be a 7-tuple.");

		/* if needed, alloc more tuples to 'out' */
		if (out.size == out.max_size)
			enlarge_ntuple_list(out);
		if (out.values == null)
			error("add_7tuple: invalid n-tuple input.");

		/* add new 7-tuple */
		out.values[out.size * out.dim + 0] = v1;
		out.values[out.size * out.dim + 1] = v2;
		out.values[out.size * out.dim + 2] = v3;
		out.values[out.size * out.dim + 3] = v4;
		out.values[out.size * out.dim + 4] = v5;
		out.values[out.size * out.dim + 5] = v6;
		out.values[out.size * out.dim + 6] = v7;

		/* update number of tuples counter */
		out.size++;
	}

	/*----------------------------------------------------------------------------*/
	/*---------------------------- 图像数据类型 -----------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	* char图像数据类型
	*
	*（x，y）处的像素值可通过以下方式访问：
	*
	* image.data [x + y * image.xsize]
	*
	*用x和y整数。
	*/
	class image_char {
		int[] data;
		int xsize, ysize;

		image_char(int xsize, int ysize) {
			this.xsize = xsize;
			this.ysize = ysize;
		}

		image_char(int xsize, int ysize, int fill_value) {
			data = new int[xsize * ysize]; /* create image */
			int N = xsize * ysize;
			int i;

			/* initialize */
			for (i = 0; i < N; i++)
				data[i] = fill_value;
			this.xsize = xsize;
			this.ysize = ysize;
		}
	}

	class image_int {
		int[] data;
		int xsize, ysize;

		image_int(int xsize, int ysize) {
			this.xsize = xsize;
			this.ysize = ysize;
		}

		image_int(int xsize, int ysize, int fill_value) {
			data = new int[xsize * ysize]; /* create image */
			int N = xsize * ysize;
			int i;

			/* 初始化 */
			for (i = 0; i < N; i++)
				data[i] = fill_value;

		}
	}

	class rect {
		double x1, y1, x2, y2; /*线段的第一个和第二个点*/
		double width; //*矩形宽度* /
		double x, y; /* 矩形的中心* */
		double theta; /* angle */
		double dx, dy; /* dx，dy）是面向矢量的线段 */
		double prec; /* 公差角度* */
		double p; /* 角度在'prec'内的点的概率 */
	};

	image_double new_image_double_ptr(int xsize, int ysize, double[] data) {

		image_double image = new image_double(xsize, ysize);

		/* check parameters */
		if (xsize == 0 || ysize == 0)
			error("new_image_double_ptr: invalid image size.");

		/* set image */
		image.data = data;

		return image;
	}

	/*----------------------------------------------------------------------------*/
	/*----------------------------- NFA computation ------------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	 * Computes the natural logarithm of the absolute value of the gamma
	 * function of x using the Lanczos approximation. See
	 * http://www.rskey.org/gamma.htm
	 * 
	 * The formula used is
	 * 
	 * @f[ \Gamma(x) = \frac{ \sum_{n=0}^{N} q_n x^n }{ \Pi_{n=0}^{N} (x+n) }
	 *     (x+5.5)^{x+0.5} e^{-(x+5.5)}
	 * @f] so
	 * @f[ \log\Gamma(x) = \log\left( \sum_{n=0}^{N} q_n x^n \right) + (x+0.5)
	 *     \log(x+5.5) - (x+5.5) - \sum_{n=0}^{N} \log(x+n)
	 * @f] and q0 = 75122.6331530, q1 = 80916.6278952, q2 = 36308.2951477, q3 =
	 *     8687.24529705, q4 = 1168.92649479, q5 = 83.8676043424, q6 =
	 *     2.50662827511.
	 */
	double log_gamma_lanczos(double x) {
		double[] q = { 75122.6331530, 80916.6278952, 36308.2951477,
				8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511 };
		double a = (x + 0.5) * Math.log(x + 5.5) - (x + 5.5);
		double b = 0.0;
		int n;

		for (n = 0; n < 7; n++) {
			a -= Math.log(x + (double) n);
			b += q[n] * Math.pow(x, (double) n);
		}
		return a + Math.log(b);
	}

	/*----------------------------------------------------------------------------*/
	 


	/**
	*计算伽玛绝对值的自然对数
	*使用Windschitl方法的x函数。 见http://www.rskey.org/gamma.htm
	*
	*使用的公式是
	*
	* @f [\ Gamma（x）= \ sqrt {\ frac {2 \ pi} {x}} \ left（\ frac {x} {e} \ sqrt {
	* x \ sinh（1 / x）+ \ frac {1} {810x ^ 6}} \ right）^ x
	* @f]如此
	* @f [\ log \ Gamma（x）= 0.5 \ log（2 \ pi）+（x-0.5）\ log（x） - x + 0.5x \ log \ left
	* x \ sinh（1 / x）+ \ frac {1} {810x ^ 6} \ right）。
	* @f]当x> 15时，这个公式是一个很好的近似值。
	*/
	double log_gamma_windschitl(double x) {
		return 0.918938533204673
				+ (x - 0.5)
				* Math.log(x)
				- x
				+ 0.5
				* x
				* Math.log(x * Math.sinh(1 / x) + 1
						/ (810.0 * Math.pow(x, 6.0)));
	}

	/*----------------------------------------------------------------------------*/
	/**
	*计算伽玛绝对值的自然对数
	*的功能。 当x> 15时使用log_gamma_windschitl（），否则使用
	* log_gamma_lanczos（）。
	*/
	double log_gamma(double x) {
		return ((x) > 15.0 ? log_gamma_windschitl(x) : log_gamma_lanczos(x));
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * 
                用于存储已计算的反转值的表的大小。
	 */
	static int TABSIZE = 100000;

	/*----------------------------------------------------------------------------*/
	/**
	 / **
*计算-log10（NFA）。
*
* NFA代表虚假报警数量：
*
* @f [\ mathrm {NFA} = NT \ cdot B（n，k，p）
* @F]
*
* - NT - 测试次数 - B（n，k，p） - 二项分布的尾部
*带参数n，k和p：
* @f [B（n，k，p）= \ sum_ {j = k} ^ n \ left（\ begin {array} {c} n \\ j \ end {array} \ right）
* p ^ {j}（1-p）^ {n-j}
* @F]
*
*值-log10（NFA）是等价的，但比NFA更直观： - -1
*对应于10个平均虚假警报 - 0对应于1个平均虚假
*警报 - 1对应于0.1平均错误警报 - 2对应于
* 0.01表示错误警报 - ...
*
*以这种方式使用，值越大，检测越好，并且a
*使用对数刻度。
* @param n
*，k，p二项参数。
* @参数logNT
*测试次数的对数
*
*计算基于伽玛函数
*以下关系：
* @f [\ left（\ begin {array} {c} n \\ k \ end {array} \ right）= \ frac {\ Gamma（n + 1）} {
* \ Gamma（k + 1）\ cdot \ Gamma（n-k + 1）}。
* @f]我们使用高效的算法来计算伽玛的对数
*功能。
*
*为了使计算速度更快，并非所有的总和都被计算出来，它的一部分
*基于对所获得的误差的限制忽略术语（an
*结果中有10％的错误被接受）。
* /
	 */
	double nfa(int n, int k, double p, double logNT) {
		
		double tolerance = 0.1; /* an error of 10% in the result is accepted */
		double log1term, term, bin_term, mult_term, bin_tail, err, p_term;
		int i;

		/* check parameters */
		if (n < 0 || k < 0 || k > n || p <= 0.0 || p >= 1.0)
			error("nfa: wrong n, k or p values.");

		/* trivial cases */
		if (n == 0 || k == 0)
			return -logNT;
		if (n == k)
			return -logNT - (double) n * Math.log10(p);

		/* probability term */
		p_term = p / (1.0 - p);

		/* compute the first term of the series */
		/*
		* binomial_tail（n，k，p）= sum_ {i = k} ^ n bincoef（n，i）* p ^ i *（1-p）^ {n-i}
        * bincoef（n，i）是二项系数。 但是bincoef（n，k）=
        * gamma（n + 1）/（gamma（k + 1）* gamma（n-k + 1））。 我们用它来计算
        *第一学期。 其实它的日志。
		 */
		log1term = log_gamma((double) n + 1.0) - log_gamma((double) k + 1.0)
				- log_gamma((double) (n - k) + 1.0) + (double) k * Math.log(p)
				+ (double) (n - k) * Math.log(1.0 - p);
		term = Math.exp(log1term);

		/* 在某些情况下不需要更多的计算 */
		if (double_equal(term, 0.0)) /* the first term is almost zero */
		{
			if ((double) k > (double) n * p) /* at begin or end of the tail? */
				return -log1term / M_LN10 - logNT; /*
													 * end: use just the first
													 * term
													 */
			else
				return -logNT; /* begin: the tail is roughly 1 */
		}

		/* compute more terms if needed */
		bin_tail = term;
		for (i = k + 1; i <= n; i++) {
			/*
			 * As term_i = bincoef(n,i) * p^i * (1-p)^(n-i) and
			 * bincoef(n,i)/bincoef(n,i-1) = n-1+1 / i, then, term_i / term_i-1
			 * = (n-i+1)/i * p/(1-p) and term_i = term_i-1 * (n-i+1)/i *
			 * p/(1-p). 1/i is stored in a table as they are computed, because
			 * divisions are expensive. p/(1-p) is computed only once and stored
			 * in 'p_term'.
			 */
			bin_term = (double) (n - i + 1)
					* (i < TABSIZE ? (inv[i] != 0.0 ? inv[i]
							: (inv[i] = 1.0 / (double) i)) : 1.0 / (double) i);

			mult_term = bin_term * p_term;
			term *= mult_term;
			bin_tail += term;
			if (bin_term < 1.0) {
				/*
				 * When bin_term<1 then mult_term_j<mult_term_i for j>i. Then,
				 * the error on the binomial tail when truncated at the i term
				 * can be bounded by a geometric series of form term_i * sum
				 * mult_term_i^j.
				 */
				err = term
						* ((1.0 - Math.pow(mult_term, (double) (n - i + 1)))
								/ (1.0 - mult_term) - 1.0);

				/*
				*最多需要容忍错误* final_result，或者：
                * tolerance * abs（-log10（bin_tail）-logNT）。 现在，这个错误
                *可以被接受在bin_tail上给出
                * tolerance * final_result除以-log10（x）的导数
                *当x = bin_tail时。 即：容忍度*
                * abs（-log10（bin_tail）-logNT）/（1 / bin_tail）最后，我们
                *如果错误小于：则截断尾部：容差*
                * abs（-log10（bin_tail）-logNT）* bin_tail
				 */
				if (err < tolerance * Math.abs(-Math.log10(bin_tail) - logNT)
						* bin_tail)
					break;
			}
		}
		return -Math.log10(bin_tail) - logNT;
	}

	/*----------------------------------------------------------------------------*/
	/*----------------------------- Gaussian filter ------------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	*计算长度为'kernel-> dim'的高斯核，标准偏差
   *'西格玛'，并且以“平均值”为中心。
   *
*例如，如果平均值= 0.5，则高斯将居中在中间
*指向值的内核 - >值[0]和内核 - >值[1]。
	 */
	void gaussian_kernel(ntuple_list kernel, double sigma, double mean) {
		double sum = 0.0;
		double val;
		int i;

		/* check parameters */
		if (kernel == null || kernel.values == null)
			error("gaussian_kernel: invalid n-tuple 'kernel'.");
		if (sigma <= 0.0)
			error("gaussian_kernel: 'sigma' must be positive.");

		/* compute Gaussian kernel */
		if (kernel.max_size < 1)
			enlarge_ntuple_list(kernel);

		kernel.size = 1;
		for (i = 0; i < kernel.dim; i++) {
			val = ((double) i - mean) / sigma;
			kernel.values[i] = Math.exp(-0.5 * val * val);
			sum += kernel.values[i];
		}

		/* normalization */
		if (sum >= 0.0)
			for (i = 0; i < kernel.dim; i++)
				kernel.values[i] /= sum;

	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Scale the input image 'in' by a factor 'scale' by Gaussian sub-sampling.
	 * 
	 * For example, scale=0.8 will give a result at 80% of the original size.
	 * 
	 * The image is convolved with a Gaussian kernel
	 * 
	 * @f[ G(x,y) = \frac{1}{2\pi\sigma^2} e^{-\frac{x^2+y^2}{2\sigma^2}}
	 * @f] before the sub-sampling to prevent aliasing.
	 * 
	 *     The standard deviation sigma given by: - sigma = sigma_scale / scale,
	 *     if scale < 1.0 - sigma = sigma_scale, if scale >= 1.0
	 * 
	 *     To be able to sub-sample at non-integer steps, some interpolation is
	 *     needed. In this implementation, the interpolation is done by the
	 *     Gaussian kernel, so both operations (filtering and sampling) are done
	 *     at the same time. The Gaussian kernel is computed centered on the
	 *     coordinates of the required sample. In this way, when applied, it
	 *     gives directly the result of convolving the image with the kernel and
	 *     interpolated to that particular position.
	 * 
	 *     A fast algorithm is done using the separability of the Gaussian
	 *     kernel. Applying the 2D Gaussian kernel is equivalent to applying
	 *     first a horizontal 1D Gaussian kernel and then a vertical 1D Gaussian
	 *     kernel (or the other way round). The reason is that
	 * @f[ G(x,y) = G(x) * G(y)
	 * @f] where
	 * @f[ G(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{x^2}{2\sigma^2}}.
	 * @f] The algorithm first applies a combined Gaussian kernel and sampling
	 *     in the x axis, and then the combined Gaussian kernel and sampling in
	 *     the y axis.
	 */
	image_double gaussian_sampler(image_double in, double scale,
			double sigma_scale) {
		image_double aux, out;
		ntuple_list kernel;
		int N, M, h, n, x, y, i;
		int xc, yc, j, double_x_size, double_y_size;
		double sigma, xx, yy, sum, prec;

		/* check parameters */
		if (in == null || in.data == null || in.xsize == 0 || in.ysize == 0)
			error("gaussian_sampler: invalid image.");
		if (scale <= 0.0)
			error("gaussian_sampler: 'scale' must be positive.");
		if (sigma_scale <= 0.0)
			error("gaussian_sampler: 'sigma_scale' must be positive.");

		/* compute new image size and get memory for images */
		if (in.xsize * scale > (double) Integer.MAX_VALUE
				|| in.ysize * scale > (double) Integer.MAX_VALUE)
			error("gaussian_sampler: the output image size exceeds the handled size.");
		N = (int) Math.ceil(in.xsize * scale);
		M = (int) Math.ceil(in.ysize * scale);
		aux = new image_double(N, in.ysize);
		out = new image_double(N, M);

		/* sigma, kernel size and memory for the kernel */
		sigma = scale < 1.0 ? sigma_scale / scale : sigma_scale;
		/*
		 * The size of the kernel is selected to guarantee that the the first
		 * discarded term is at least 10^prec times smaller than the central
		 * value. For that, h should be larger than x, with e^(-x^2/2sigma^2) =
		 * 1/10^prec. Then, x = sigma * sqrt( 2 * prec * ln(10) ).
		 */
		prec = 3.0;
		h = (int) Math.ceil(sigma * Math.sqrt(2.0 * prec * Math.log(10.0)));
		n = 1 + 2 * h; /* kernel size */
		kernel = new ntuple_list(n);

		/* auxiliary double image size variables */
		double_x_size = (int) (2 * in.xsize);
		double_y_size = (int) (2 * in.ysize);

		/* First subsampling: x axis */
		for (x = 0; x < aux.xsize; x++) {
			/*
			 * x is the coordinate in the new image. xx is the corresponding
			 * x-value in the original size image. xc is the integer value, the
			 * pixel coordinate of xx.
			 */
			xx = (double) x / scale;
			/*
			 * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
			 * pixel with xc=0 get the values of xx from -0.5 to 0.5
			 */
			xc = (int) Math.floor(xx + 0.5);
			gaussian_kernel(kernel, sigma, (double) h + xx - (double) xc);
			/*
			 * the kernel must be computed for each x because the fine offset
			 * xx-xc is different in each case
			 */

			for (y = 0; y < aux.ysize; y++) {
				sum = 0.0;
				for (i = 0; i < kernel.dim; i++) {
					j = xc - h + i;

					/* symmetry boundary condition */
					while (j < 0)
						j += double_x_size;
					while (j >= double_x_size)
						j -= double_x_size;
					if (j >= (int) in.xsize)
						j = double_x_size - 1 - j;

					sum += in.data[j + y * in.xsize] * kernel.values[i];
				}
				aux.data[x + y * aux.xsize] = sum;
			}
		}

		/* Second subsampling: y axis */
		for (y = 0; y < out.ysize; y++) {
			/*
			 * y is the coordinate in the new image. yy is the corresponding
			 * x-value in the original size image. yc is the integer value, the
			 * pixel coordinate of xx.
			 */
			yy = (double) y / scale;
			/*
			 * coordinate (0.0,0.0) is in the center of pixel (0,0), so the
			 * pixel with yc=0 get the values of yy from -0.5 to 0.5
			 */
			yc = (int) Math.floor(yy + 0.5);
			gaussian_kernel(kernel, sigma, (double) h + yy - (double) yc);
			/*
			 * the kernel must be computed for each y because the fine offset
			 * yy-yc is different in each case
			 */

			for (x = 0; x < out.xsize; x++) {
				sum = 0.0;
				for (i = 0; i < kernel.dim; i++) {
					j = yc - h + i;

					/* symmetry boundary condition */
					while (j < 0)
						j += double_y_size;
					while (j >= double_y_size)
						j -= double_y_size;
					if (j >= (int) in.ysize)
						j = double_y_size - 1 - j;

					sum += aux.data[x + j * aux.xsize] * kernel.values[i];
				}
				out.data[x + y * out.xsize] = sum;
			}
		}

		return out;
	}

	/*----------------------------------------------------------------------------*/
	/*--------------------------------- Gradient ---------------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	 * Computes the direction of the level line of 'in' at each point.
	 * 
	 * The result is: - an image_double with the angle at each pixel, or NOTDEF
	 * if not defined. - the image_double 'modgrad' (a pointer is passed as
	 * argument) with the gradient magnitude at each point. - a list of pixels
	 * 'list_p' roughly ordered by decreasing gradient magnitude. (The order is
	 * made by classifying points into bins by gradient magnitude. The
	 * parameters 'n_bins' and 'max_grad' specify the number of bins and the
	 * gradient modulus at the highest bin. The pixels in the list would be in
	 * decreasing gradient magnitude, up to a precision of the size of the
	 * bins.) - a pointer 'mem_p' to the memory used by 'list_p' to be able to
	 * free the memory when it is not used anymore.
	 */

	image_double modgrad;
	coorlist[] mem_p;
	coorlist list_p;

	image_double ll_angle(image_double in, double threshold,
	/* coorlist [] list_p, *//* void ** mem_p, */
	/* image_double modgrad, */int n_bins) {
		image_double g;
		int n, p, x, y, adr, i;
		double com1, com2, gx, gy, norm, norm2;
		/*
		 * the rest of the variables are used for pseudo-ordering the gradient
		 * magnitude values
		 */
		int list_count = 0;
		coorlist[] list;

		coorlist[] range_l_s; /* array of pointers to start of bin list */
		coorlist[] range_l_e; /* array of pointers to end of bin list */
		coorlist start;
		coorlist end;
		double max_grad = 0.0;

		/* check parameters */
		if (in == null || in.data == null || in.xsize == 0 || in.ysize == 0)
			error("ll_angle: invalid image.");
		if (threshold < 0.0)
			error("ll_angle: 'threshold' must be positive.");
		if (list_p == null) {
			error("ll_angle: null pointer 'list_p'.");
			// list_p = new coorlist();
		}
		// if (mem_p == null)
		// error("ll_angle: null pointer 'mem_p'.");
		// if (modgrad == null)
		// error("ll_angle: null pointer 'modgrad'.");
		if (n_bins == 0)
			error("ll_angle: 'n_bins' must be positive.");

		/* image size shortcuts */
		n = in.ysize;
		p = in.xsize;

		list = new coorlist[n * p];
		for (int z = 0; z < n * p; z++) {
			list[z] = new coorlist();
		}

		mem_p = list;

		/* allocate output image */
		g = new image_double(in.xsize, in.ysize);

		/* get memory for the image of gradient modulus */
		modgrad = new image_double(in.xsize, in.ysize);

		/* get memory for "ordered" list of pixels */
		// list = new coorlist[n * p];// (struct coorlist *) calloc( (size_t)
		// (n*p), sizeof(struct coorlist) );

		range_l_s = new coorlist[n_bins];
		range_l_e = new coorlist[n_bins];

		if (list == null || range_l_s == null || range_l_e == null)
			error("not enough memory.");

		for (i = 0; i < n_bins; i++) {
			range_l_s[i] = range_l_e[i] = null;
		}

		/* 'undefined' on the down and right boundaries */
		for (x = 0; x < p; x++)
			g.data[(n - 1) * p + x] = NOTDEF;

		for (y = 0; y < n; y++)
			g.data[p * y + p - 1] = NOTDEF;

		/* compute gradient on the remaining pixels */
		for (x = 0; x < p - 1; x++)
			for (y = 0; y < n - 1; y++) {
				adr = y * p + x;

				/*
				 * Norm 2 computation using 2x2 pixel window: A B C D and com1 =
				 * D-A, com2 = B-C. Then gx = B+D - (A+C) horizontal difference
				 * gy = C+D - (A+B) vertical difference com1 and com2 are just
				 * to avoid 2 additions.
				 */
				com1 = in.data[adr + p + 1] - in.data[adr];
				com2 = in.data[adr + 1] - in.data[adr + p];

				gx = com1 + com2; /* gradient x component */
				gy = com1 - com2; /* gradient y component */
				norm2 = gx * gx + gy * gy;
				norm = Math.sqrt(norm2 / 4.0); /* gradient norm */

				modgrad.data[adr] = norm; /* store gradient norm */

				// System.out.println("norm " + norm + " threshold " +
				// threshold);

				if (norm <= threshold) /* norm too small, gradient no defined */
					g.data[adr] = NOTDEF; /* gradient angle not defined */
				else {
					// System.out.println("gradient angle --------------");
					/* gradient angle computation */
					g.data[adr] = Math.atan2(gx, -gy);

					/* look for the maximum of the gradient */
					if (norm > max_grad)
						max_grad = norm;
				}
			}

		/* compute histogram of gradient values */
		for (x = 0; x < p - 1; x++)
			for (y = 0; y < n - 1; y++) {
				norm = modgrad.data[y * p + x];

				/* store the point in the right bin according to its norm */
				i = (int) (norm * (double) n_bins / max_grad);
				if (i >= n_bins)
					i = n_bins - 1;
				if (range_l_e[i] == null) {
					// System.out.println("here1 "+list[list_count]);

					range_l_s[i] = range_l_e[i] = list[list_count++];
				} else {
					// System.out.println("here2");
					range_l_e[i].next = list[list_count];
					range_l_e[i] = list[list_count++];
				}
				// System.out.println(i + " "+range_l_e[i]);
				range_l_e[i].x = (int) x;
				range_l_e[i].y = (int) y;
				range_l_e[i].next = null;
			}

		/*
		 * Make the list of pixels (almost) ordered by norm value. It starts by
		 * the larger bin, so the list starts by the pixels with the highest
		 * gradient value. Pixels would be ordered by norm value, up to a
		 * precision given by max_grad/n_bins.
		 */
		for (i = n_bins - 1; i > 0 && range_l_s[i] == null; i--) {

		}

	//	System.out.println("i val " + i);

		start = range_l_s[i];
		end = range_l_e[i];
		if (start != null) {
			//System.out.println("start not null");
			while (i > 0) {
				--i;
				if (range_l_s[i] != null) {
					// System.out.println("range not null");

					end.next = range_l_s[i];
					end = range_l_e[i];
				}
			}
		}
		list_p = start;

		return g;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Absolute value angle difference.
	 */
	double angle_diff(double a, double b) {
		a -= b;
		while (a <= -M_PI)
			a += M_2__PI;
		while (a > M_PI)
			a -= M_2__PI;
		if (a < 0.0)
			a = -a;
		return a;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Signed angle difference.
	 */
	double angle_diff_signed(double a, double b) {
		a -= b;
		while (a <= -M_PI)
			a += M_2__PI;
		while (a > M_PI)
			a -= M_2__PI;
		return a;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Is point (x,y) aligned to angle theta, up to precision 'prec'?
	 */
	boolean isaligned(int x, int y, image_double angles, double theta,
			double prec) {
		double a;

		/* check parameters */
		if (angles == null || angles.data == null)
			error("isaligned: invalid image 'angles'.");
		if (x < 0 || y < 0 || x >= (int) angles.xsize
				|| y >= (int) angles.ysize)
			error("isaligned: (x,y) out of the image.");
		if (prec < 0.0)
			error("isaligned: 'prec' must be positive.");

		/* angle at pixel (x,y) */
		a = angles.data[x + y * angles.xsize];

		/*
		 * pixels whose level-line angle is not defined are considered as
		 * NON-aligned
		 */
		if (a == NOTDEF)
			return false; /*
						 * there is no need to call the function 'double_equal'
						 * here because there is no risk of problems related to
						 * the comparison doubles, we are only interested in the
						 * exact NOTDEF value
						 */

		/* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
		theta -= a;
		if (theta < 0.0)
			theta = -theta;
		if (theta > M_3_2_PI) {
			theta -= M_2__PI;
			if (theta < 0.0)
				theta = -theta;
		}

		return theta <= prec;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Build a region of pixels that share the same angle, up to a tolerance
	 * 'prec', starting at point (x,y).
	 */

	public int reg_size;
	public double reg_angle;

	void region_grow(int x, int y, image_double angles, Point[] reg,
	/* int * reg_size, double * reg_angle, */image_char used, double prec) {
		double sumdx, sumdy;
		int xx, yy, i;

		/* check parameters */
		if (x < 0 || y < 0 || x >= (int) angles.xsize
				|| y >= (int) angles.ysize)
			error("region_grow: (x,y) out of the image.");
		if (angles == null || angles.data == null)
			error("region_grow: invalid image 'angles'.");
		if (reg == null)
			error("region_grow: invalid 'reg'.");
		// if( reg_size == null )
		// error("region_grow: invalid pointer 'reg_size'.");
		// if( reg_angle == NULL )
		// error("region_grow: invalid pointer 'reg_angle'.");
		// if( used == NULL || used->data == NULL )
		// error("region_grow: invalid image 'used'.");

		/* first point of the region */
		reg_size = 1;
		reg[0].x = x;
		reg[0].y = y;
		reg_angle = angles.data[x + y * angles.xsize]; /* region's angle */
		sumdx = Math.cos(reg_angle);
		sumdy = Math.sin(reg_angle);
		used.data[x + y * used.xsize] = USED;

		/* try neighbors as new region points */
		for (i = 0; i < reg_size; i++)
			for (xx = reg[i].x - 1; xx <= reg[i].x + 1; xx++)
				for (yy = reg[i].y - 1; yy <= reg[i].y + 1; yy++) {
					
					//System.out.println("here2____"+(xx >= 0 && yy >= 0 && xx < (int) used.xsize
					//		&& yy < (int) used.ysize)/*
					//		&& used.data[xx + yy * used.xsize] != USED*/+"__________");
					
					//System.out.println("::::::::: "+xx+","+yy+"  "+used.xsize+","+used.ysize);
					
					if (xx >= 0 && yy >= 0 && xx < (int) used.xsize
							&& yy < (int) used.ysize
							&& used.data[xx + yy * used.xsize] != USED
							&& isaligned(xx, yy, angles, reg_angle, prec)) {
						/* add point */
						used.data[xx + yy * used.xsize] = USED;
						reg[reg_size].x = xx;
						reg[reg_size].y = yy;
						++(reg_size);

						/* update region's angle */
						sumdx += Math.cos(angles.data[xx + yy * angles.xsize]);
						sumdy += Math.sin(angles.data[xx + yy * angles.xsize]);
						reg_angle = Math.atan2(sumdy, sumdx);
					}
				}
		//System.out.println(">>>regsize " + reg_size);

	}

	/*----------------------------------------------------------------------------*/
	/*---------------------------------- Regions ---------------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	 * Compute region's angle as the principal inertia axis of the region.
	 * 
	 * The following is the region inertia matrix A:
	 * 
	 * @f[
	 * 
	 *     A = \left(\begin{array}{cc} Ixx & Ixy \\ Ixy & Iyy \\
	 *     \end{array}\right)
	 * @f] where
	 * 
	 *     Ixx = sum_i G(i).(y_i - cx)^2
	 * 
	 *     Iyy = sum_i G(i).(x_i - cy)^2
	 * 
	 *     Ixy = - sum_i G(i).(x_i - cx).(y_i - cy)
	 * 
	 *     and - G(i) is the gradient norm at pixel i, used as pixel's weight. -
	 *     x_i and y_i are the coordinates of pixel i. - cx and cy are the
	 *     coordinates of the center of th region.
	 * 
	 *     lambda1 and lambda2 are the eigenvalues of matrix A, with lambda1 >=
	 *     lambda2. They are found by solving the characteristic polynomial:
	 * 
	 *     det( lambda I - A) = 0
	 * 
	 *     that gives:
	 * 
	 *     lambda1 = ( Ixx + Iyy + sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2
	 * 
	 *     lambda2 = ( Ixx + Iyy - sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2
	 * 
	 *     To get the line segment direction we want to get the angle the
	 *     eigenvector associated to the smallest eigenvalue. We have to solve
	 *     for a,b in:
	 * 
	 *     a.Ixx + b.Ixy = a.lambda2
	 * 
	 *     a.Ixy + b.Iyy = b.lambda2
	 * 
	 *     We want the angle theta = atan(b/a). It can be computed with any of
	 *     the two equations:
	 * 
	 *     theta = atan( (lambda2-Ixx) / Ixy )
	 * 
	 *     or
	 * 
	 *     theta = atan( Ixy / (lambda2-Iyy) )
	 * 
	 *     When |Ixx| > |Iyy| we use the first, otherwise the second (just to
	 *     get better numeric precision).
	 */
	double get_theta(Point[] reg, int reg_size, double x, double y,
			image_double modgrad, double reg_angle, double prec) {
		double lambda, theta, weight;
		double Ixx = 0.0;
		double Iyy = 0.0;
		double Ixy = 0.0;
		int i;

		/* check parameters */
		if (reg == null)
			error("get_theta: invalid region.");
		if (reg_size <= 1)
			error("get_theta: region size <= 1.");
		if (modgrad == null || modgrad.data == null)
			error("get_theta: invalid 'modgrad'.");
		if (prec < 0.0)
			error("get_theta: 'prec' must be positive.");

		/* compute inertia matrix */
		for (i = 0; i < reg_size; i++) {
			weight = modgrad.data[reg[i].x + reg[i].y * modgrad.xsize];
			Ixx += ((double) reg[i].y - y) * ((double) reg[i].y - y) * weight;
			Iyy += ((double) reg[i].x - x) * ((double) reg[i].x - x) * weight;
			Ixy -= ((double) reg[i].x - x) * ((double) reg[i].y - y) * weight;
		}
		if (double_equal(Ixx, 0.0) && double_equal(Iyy, 0.0)
				&& double_equal(Ixy, 0.0))
			error("get_theta: null inertia matrix.");

		/* compute smallest eigenvalue */
		lambda = 0.5 * (Ixx + Iyy - Math.sqrt((Ixx - Iyy) * (Ixx - Iyy) + 4.0
				* Ixy * Ixy));

		/* compute angle */
		theta = Math.abs(Ixx) > Math.abs(Iyy) ? Math.atan2(lambda - Ixx, Ixy)
				: Math.atan2(Ixy, lambda - Iyy);

		/*
		 * The previous procedure doesn't cares about orientation, so it could
		 * be wrong by 180 degrees. Here is corrected if necessary.
		 */
		if (angle_diff(theta, reg_angle) > prec)
			theta += M_PI;

		return theta;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Reduce the region size, by elimination the points far from the starting
	 * point, until that leads to rectangle with the right density of region
	 * points or to discard the region if too small.
	 */
	boolean reduce_region_radius(Point[] reg, /* int * reg_size, */
			image_double modgrad, /*double reg_angle,*/ double prec, double p,
			rect rec, image_char used, image_double angles, double density_th) {
		double density, rad1, rad2, rad, xc, yc;
		int i;

		/* check parameters */
		if (reg == null)
			error("reduce_region_radius: invalid pointer 'reg'.");
		// if( reg_size == null )
		// error("reduce_region_radius: invalid pointer 'reg_size'.");
		if (prec < 0.0)
			error("reduce_region_radius: 'prec' must be positive.");
		if (rec == null)
			error("reduce_region_radius: invalid pointer 'rec'.");
		if (used == null || used.data == null)
			error("reduce_region_radius: invalid image 'used'.");
		if (angles == null || angles.data == null)
			error("reduce_region_radius: invalid image 'angles'.");

		/* compute region points density */
		density = (double) reg_size
				/ (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

		/* if the density criterion is satisfied there is nothing to do */
		if (density >= density_th)
			return true;

		/* compute region's radius */
		xc = (double) reg[0].x;
		yc = (double) reg[0].y;
		rad1 = dist(xc, yc, rec.x1, rec.y1);
		rad2 = dist(xc, yc, rec.x2, rec.y2);
		rad = rad1 > rad2 ? rad1 : rad2;

		/* while the density criterion is not satisfied, remove farther pixels */
		while (density < density_th) {
			rad *= 0.75; /* reduce region's radius to 75% of its value */

			/* remove points from the region and update 'used' map */
			for (i = 0; i < reg_size; i++)
				if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) > rad) {
					/* point not kept, mark it as NOTUSED */
					used.data[reg[i].x + reg[i].y * used.xsize] = NOTUSED;
					/* remove point from the region */
					reg[i].x = reg[reg_size - 1].x; /*
													 * if i==*reg_size-1 copy
													 * itself
													 */
					reg[i].y = reg[reg_size - 1].y;
					--(reg_size);
					--i; /* to avoid skipping one point */
				}

			/*
			 * reject if the region is too small. 2 is the minimal region size
			 * for 'region2rect' to work.
			 */
			if (reg_size < 2)
				return false;

			/* re-compute rectangle */
			region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);

			/* re-compute region points density */
			density = (double) reg_size
					/ (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);
		}

		/* if this point is reached, the density criterion is satisfied */
		return true;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Computes a rectangle that covers a region of points.
	 */
	void region2rect(Point[] reg, int reg_size, image_double modgrad,
			/*double reg_angle,*/ double prec, double p, rect rec) {
		double x, y, dx, dy, l, w, theta, weight, sum, l_min, l_max, w_min, w_max;
		int i;

		/* check parameters */
		if (reg == null)
			error("region2rect: invalid region.");
		if (reg_size <= 1)
			error("region2rect: region size <= 1.");
		if (modgrad == null || modgrad.data == null)
			error("region2rect: invalid image 'modgrad'.");
		if (rec == null)
			error("region2rect: invalid 'rec'.");

		/*
		 * center of the region:
		 * 
		 * It is computed as the weighted sum of the coordinates of all the
		 * pixels in the region. The norm of the gradient is used as the weight
		 * of a pixel. The sum is as follows: cx = \sum_i G(i).x_i cy = \sum_i
		 * G(i).y_i where G(i) is the norm of the gradient of pixel i and
		 * x_i,y_i are its coordinates.
		 */
		x = y = sum = 0.0;
		for (i = 0; i < reg_size; i++) {
			weight = modgrad.data[reg[i].x + reg[i].y * modgrad.xsize];
			x += (double) reg[i].x * weight;
			y += (double) reg[i].y * weight;
			sum += weight;
		}
		if (sum <= 0.0)
			error("region2rect: weights sum equal to zero.");
		x /= sum;
		y /= sum;

		/* theta */
		theta = get_theta(reg, reg_size, x, y, modgrad, reg_angle, prec);

		/*
		 * length and width:
		 * 
		 * 'l' and 'w' are computed as the distance from the center of the
		 * region to pixel i, projected along the rectangle axis (dx,dy) and to
		 * the orthogonal axis (-dy,dx), respectively.
		 * 
		 * The length of the rectangle goes from l_min to l_max, where l_min and
		 * l_max are the minimum and maximum values of l in the region.
		 * Analogously, the width is selected from w_min to w_max, where w_min
		 * and w_max are the minimum and maximum of w for the pixels in the
		 * region.
		 */
		dx = Math.cos(theta);
		dy = Math.sin(theta);
		l_min = l_max = w_min = w_max = 0.0;
		for (i = 0; i < reg_size; i++) {
			l = ((double) reg[i].x - x) * dx + ((double) reg[i].y - y) * dy;
			w = -((double) reg[i].x - x) * dy + ((double) reg[i].y - y) * dx;

			if (l > l_max)
				l_max = l;
			if (l < l_min)
				l_min = l;
			if (w > w_max)
				w_max = w;
			if (w < w_min)
				w_min = w;
		}

		/* store values */
		rec.x1 = x + l_min * dx;
		rec.y1 = y + l_min * dy;
		rec.x2 = x + l_max * dx;
		rec.y2 = y + l_max * dy;
		rec.width = w_max - w_min;
		rec.x = x;
		rec.y = y;
		rec.theta = theta;
		rec.dx = dx;
		rec.dy = dy;
		rec.prec = prec;
		rec.p = p;

		/*
		 * we impose a minimal width of one pixel
		 * 
		 * A sharp horizontal or vertical step would produce a perfectly
		 * horizontal or vertical region. The width computed would be zero. But
		 * that corresponds to a one pixels width transition in the image.
		 */
		if (rec.width < 1.0)
			rec.width = 1.0;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Refine a rectangle.
	 * 
	 * For that, an estimation of the angle tolerance is performed by the
	 * standard deviation of the angle at points near the region's starting
	 * point. Then, a new region is grown starting from the same point, but
	 * using the estimated angle tolerance. If this fails to produce a rectangle
	 * with the right density of region points, 'reduce_region_radius' is called
	 * to try to satisfy this condition.
	 */
	boolean refine(Point[] reg, /*int reg_size,*/ image_double modgrad,
			/*double reg_angle,*/ double prec, double p, rect rec, image_char used,
			image_double angles, double density_th) {
		double angle, ang_d, mean_angle, tau, density, xc, yc, ang_c, sum, s_sum;
		int i, n;

		//this.reg_size = reg_size;
		//this.reg_angle = reg_angle;

		/* check parameters */
		if (reg == null)
			error("refine: invalid pointer 'reg'.");
		// if( reg_size == null ) error("refine: invalid pointer 'reg_size'.");
		if (prec < 0.0)
			error("refine: 'prec' must be positive.");
		if (rec == null)
			error("refine: invalid pointer 'rec'.");
		if (used == null || used.data == null)
			error("refine: invalid image 'used'.");
		if (angles == null || angles.data == null)
			error("refine: invalid image 'angles'.");

		/* compute region points density */
		density = (double) reg_size
				/ (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

		/* if the density criterion is satisfied there is nothing to do */
		if (density >= density_th)
			return true;

		/*------ First try: reduce angle tolerance ------*/

		/* compute the new mean angle and tolerance */
		xc = (double) reg[0].x;
		yc = (double) reg[0].y;
		ang_c = angles.data[reg[0].x + reg[0].y * angles.xsize];
		sum = s_sum = 0.0;
		n = 0;
		for (i = 0; i < this.reg_size; i++) {
			used.data[reg[i].x + reg[i].y * used.xsize] = NOTUSED;
			if (dist(xc, yc, (double) reg[i].x, (double) reg[i].y) < rec.width) {
				angle = angles.data[reg[i].x + reg[i].y * angles.xsize];
				ang_d = angle_diff_signed(angle, ang_c);
				sum += ang_d;
				s_sum += ang_d * ang_d;
				++n;
			}
		}
		mean_angle = sum / (double) n;
		tau = 2.0 * Math.sqrt((s_sum - 2.0 * mean_angle * sum) / (double) n
				+ mean_angle * mean_angle); /* 2 * standard deviation */

		/*
		 * find a new region from the same starting point and new angle
		 * tolerance
		 */
		region_grow(reg[0].x, reg[0].y, angles, reg /* ,reg_size,reg_angle, */,
				used, tau);
		/* if the region is too small, reject */
		if (reg_size < 2)
			return false;

		/* re-compute rectangle */
		region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);

		/* re-compute region points density */
		density = (double) reg_size
				/ (dist(rec.x1, rec.y1, rec.x2, rec.y2) * rec.width);

		/*------ Second try: reduce region radius ------*/
		if (density < density_th)
			return reduce_region_radius(reg, /* reg_size, */modgrad, /*reg_angle,*/
					prec, p, rec, used, angles, density_th);

		/* if this point is reached, the density criterion is satisfied */
		return true;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Interpolate y value corresponding to 'x' value given, in the line 'x1,y1'
	 * to 'x2,y2'; if 'x1=x2' return the smaller of 'y1' and 'y2'.
	 * 
	 * The following restrictions are required: - x1 <= x2 - x1 <= x - x <= x2
	 */
	double inter_low(double x, double x1, double y1, double x2, double y2) {
		/* check parameters */
		if (x1 > x2 || x < x1 || x > x2)
			error("inter_low: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");

		/* interpolation */
		if (double_equal(x1, x2) && y1 < y2)
			return y1;
		if (double_equal(x1, x2) && y1 > y2)
			return y2;
		return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Interpolate y value corresponding to 'x' value given, in the line 'x1,y1'
	 * to 'x2,y2'; if 'x1=x2' return the larger of 'y1' and 'y2'.
	 * 
	 * The following restrictions are required: - x1 <= x2 - x1 <= x - x <= x2
	 */
	double inter_hi(double x, double x1, double y1, double x2, double y2) {
		/* check parameters */
		if (x1 > x2 || x < x1 || x > x2)
			error("inter_hi: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");

		/* interpolation */
		if (double_equal(x1, x2) && y1 < y2)
			return y2;
		if (double_equal(x1, x2) && y1 > y2)
			return y1;
		return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Free memory used by a rectangle iterator.
	 */
	void ri_del(rect_itr iter) {
		if (iter == null)
			error("ri_del: NULL iterator.");
		// free( (void *) iter );
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Check if the iterator finished the full iteration.
	 * 
	 * See details in \ref rect_iter
	 */
	boolean ri_end(rect_itr i) {
		/* check input */
		if (i == null)
			error("ri_end: NULL iterator.");

		/*
		 * if the current x value is larger than the largest x value in the
		 * rectangle (vx[2]), we know the full exploration of the rectangle is
		 * finished.
		 */
		return (double) (i.x) > i.vx[2];
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Increment a rectangle iterator.
	 * 
	 * See details in \ref rect_iter
	 */
	void ri_inc(rect_itr i) {
		/* check input */
		if (i == null)
			error("ri_inc: NULL iterator.");

		/*
		 * if not at end of exploration, increase y value for next pixel in the
		 * 'column'
		 */
		if (!ri_end(i))
			i.y++;

		/*
		 * if the end of the current 'column' is reached, and it is not the end
		 * of exploration, advance to the next 'column'
		 */
		while ((double) (i.y) > i.ye && !ri_end(i)) {
			/* increase x, next 'column' */
			i.x++;

			/* if end of exploration, return */
			if (ri_end(i))
				return;

			/*
			 * update lower y limit (start) for the new 'column'.
			 * 
			 * We need to interpolate the y value that corresponds to the lower
			 * side of the rectangle. The first thing is to decide if the
			 * corresponding side is
			 * 
			 * vx[0],vy[0] to vx[3],vy[3] or vx[3],vy[3] to vx[2],vy[2]
			 * 
			 * Then, the side is interpolated for the x value of the 'column'.
			 * But, if the side is vertical (as it could happen if the rectangle
			 * is vertical and we are dealing with the first or last 'columns')
			 * then we pick the lower value of the side by using 'inter_low'.
			 */
			if ((double) i.x < i.vx[3])
				i.ys = inter_low((double) i.x, i.vx[0], i.vy[0], i.vx[3],
						i.vy[3]);
			else
				i.ys = inter_low((double) i.x, i.vx[3], i.vy[3], i.vx[2],
						i.vy[2]);

			/*
			 * update upper y limit (end) for the new 'column'.
			 * 
			 * We need to interpolate the y value that corresponds to the upper
			 * side of the rectangle. The first thing is to decide if the
			 * corresponding side is
			 * 
			 * vx[0],vy[0] to vx[1],vy[1] or vx[1],vy[1] to vx[2],vy[2]
			 * 
			 * Then, the side is interpolated for the x value of the 'column'.
			 * But, if the side is vertical (as it could happen if the rectangle
			 * is vertical and we are dealing with the first or last 'columns')
			 * then we pick the lower value of the side by using 'inter_low'.
			 */
			if ((double) i.x < i.vx[1])
				i.ye = inter_hi((double) i.x, i.vx[0], i.vy[0], i.vx[1],
						i.vy[1]);
			else
				i.ye = inter_hi((double) i.x, i.vx[1], i.vy[1], i.vx[2],
						i.vy[2]);

			/* new y */
			i.y = (int) Math.ceil(i.ys);
		}
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Create and initialize a rectangle iterator.
	 * 
	 * See details in \ref rect_iter
	 */
	rect_itr ri_ini(rect r) {
		double[] vx = new double[4];
		double[] vy = new double[4];
		int n, offset;
		rect_itr i;

		/* check parameters */
		if (r == null)
			error("ri_ini: invalid rectangle.");

		/* get memory */
		i = new rect_itr();
		if (i == null)
			error("ri_ini: Not enough memory.");

		/*
		 * build list of rectangle corners ordered in a circular way around the
		 * rectangle
		 */
		vx[0] = r.x1 - r.dy * r.width / 2.0;
		vy[0] = r.y1 + r.dx * r.width / 2.0;
		vx[1] = r.x2 - r.dy * r.width / 2.0;
		vy[1] = r.y2 + r.dx * r.width / 2.0;
		vx[2] = r.x2 + r.dy * r.width / 2.0;
		vy[2] = r.y2 - r.dx * r.width / 2.0;
		vx[3] = r.x1 + r.dy * r.width / 2.0;
		vy[3] = r.y1 - r.dx * r.width / 2.0;

		/*
		 * compute rotation of index of corners needed so that the first point
		 * has the smaller x.
		 * 
		 * if one side is vertical, thus two corners have the same smaller x
		 * value, the one with the largest y value is selected as the first.
		 */
		if (r.x1 < r.x2 && r.y1 <= r.y2)
			offset = 0;
		else if (r.x1 >= r.x2 && r.y1 < r.y2)
			offset = 1;
		else if (r.x1 > r.x2 && r.y1 >= r.y2)
			offset = 2;
		else
			offset = 3;

		/* apply rotation of index. */
		for (n = 0; n < 4; n++) {
			i.vx[n] = vx[(offset + n) % 4];
			i.vy[n] = vy[(offset + n) % 4];
		}

		/*
		 * Set an initial condition.
		 * 
		 * The values are set to values that will cause 'ri_inc' (that will be
		 * called immediately) to initialize correctly the first 'column' and
		 * compute the limits 'ys' and 'ye'.
		 * 
		 * 'y' is set to the integer value of vy[0], the starting corner.
		 * 
		 * 'ys' and 'ye' are set to very small values, so 'ri_inc' will notice
		 * that it needs to start a new 'column'.
		 * 
		 * The smallest integer coordinate inside of the rectangle is
		 * 'ceil(vx[0])'. The current 'x' value is set to that value minus one,
		 * so 'ri_inc' (that will increase x by one) will advance to the first
		 * 'column'.
		 */
		i.x = (int) Math.ceil(i.vx[0]) - 1;
		i.y = (int) Math.ceil(i.vy[0]);
		i.ys = i.ye = -Double.MAX_VALUE;

		/* advance to the first pixel */
		ri_inc(i);

		return i;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Compute a rectangle's NFA value.
	 */
	double rect_nfa(rect rec, image_double angles, double logNT) {
		rect_itr i;
		int pts = 0;
		int alg = 0;

		/* check parameters */
		if (rec == null)
			error("rect_nfa: invalid rectangle.");
		if (angles == null)
			error("rect_nfa: invalid 'angles'.");

		/* compute the total number of pixels and of aligned points in 'rec' */
		for (i = ri_ini(rec); !ri_end(i); ri_inc(i))
			/* rectangle iterator */
			if (i.x >= 0 && i.y >= 0 && i.x < (int) angles.xsize
					&& i.y < (int) angles.ysize) {
				++pts; /* total number of pixels counter */
				if (isaligned(i.x, i.y, angles, rec.theta, rec.prec))
					++alg; /* aligned points counter */
			}
		// ri_del(i); /* delete iterator */

		return nfa(pts, alg, rec.p, logNT); /* compute NFA value */
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Copy one rectangle structure to another.
	 */
	void rect_copy(rect in, rect out) {
		/* check parameters */
		if (in == null || out == null)
			error("rect_copy: invalid 'in' or 'out'.");

		/* copy values */
		out.x1 = in.x1;
		out.y1 = in.y1;
		out.x2 = in.x2;
		out.y2 = in.y2;
		out.width = in.width;
		out.x = in.x;
		out.y = in.y;
		out.theta = in.theta;
		out.dx = in.dx;
		out.dy = in.dy;
		out.prec = in.prec;
		out.p = in.p;
	}

	/*----------------------------------------------------------------------------*/
	/**
	 * Try some rectangles variations to improve NFA value. Only if the
	 * rectangle is not meaningful (i.e., log_nfa <= log_eps).
	 */
	double rect_improve(rect rec, image_double angles, double logNT,
			double log_eps) {
		rect r = new rect();
		double log_nfa, log_nfa_new;
		double delta = 0.5;
		double delta_2 = delta / 2.0;
		int n;

		log_nfa = rect_nfa(rec, angles, logNT);

		if (log_nfa > log_eps)
			return log_nfa;

		/* try finer precisions */
		rect_copy(rec, r);
		for (n = 0; n < 5; n++) {
			r.p /= 2.0;
			r.prec = r.p * M_PI;
			log_nfa_new = rect_nfa(r, angles, logNT);
			if (log_nfa_new > log_nfa) {
				log_nfa = log_nfa_new;
				rect_copy(r, rec);
			}
		}

		if (log_nfa > log_eps)
			return log_nfa;

		/* try to reduce width */
		rect_copy(rec, r);
		for (n = 0; n < 5; n++) {
			if ((r.width - delta) >= 0.5) {
				r.width -= delta;
				log_nfa_new = rect_nfa(r, angles, logNT);
				if (log_nfa_new > log_nfa) {
					rect_copy(r, rec);
					log_nfa = log_nfa_new;
				}
			}
		}

		if (log_nfa > log_eps)
			return log_nfa;

		/* try to reduce one side of the rectangle */
		rect_copy(rec, r);
		for (n = 0; n < 5; n++) {
			if ((r.width - delta) >= 0.5) {
				r.x1 += -r.dy * delta_2;
				r.y1 += r.dx * delta_2;
				r.x2 += -r.dy * delta_2;
				r.y2 += r.dx * delta_2;
				r.width -= delta;
				log_nfa_new = rect_nfa(r, angles, logNT);
				if (log_nfa_new > log_nfa) {
					rect_copy(r, rec);
					log_nfa = log_nfa_new;
				}
			}
		}

		if (log_nfa > log_eps)
			return log_nfa;

		/* try to reduce the other side of the rectangle */
		rect_copy(rec, r);
		for (n = 0; n < 5; n++) {
			if ((r.width - delta) >= 0.5) {
				r.x1 -= -r.dy * delta_2;
				r.y1 -= r.dx * delta_2;
				r.x2 -= -r.dy * delta_2;
				r.y2 -= r.dx * delta_2;
				r.width -= delta;
				log_nfa_new = rect_nfa(r, angles, logNT);
				if (log_nfa_new > log_nfa) {
					rect_copy(r, rec);
					log_nfa = log_nfa_new;
				}
			}
		}

		if (log_nfa > log_eps)
			return log_nfa;

		/* try even finer precisions */
		rect_copy(rec, r);
		for (n = 0; n < 5; n++) {
			r.p /= 2.0;
			r.prec = r.p * M_PI;
			log_nfa_new = rect_nfa(r, angles, logNT);
			if (log_nfa_new > log_nfa) {
				log_nfa = log_nfa_new;
				rect_copy(r, rec);
			}
		}

		return log_nfa;
	}

	/*----------------------------------------------------------------------------*/
	/*-------------------------- Line Segment Detector ---------------------------*/
	/*----------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------*/
	/**
	 * LSD full interface.
	 */
	double[] LineSegmentDetection(double[] img, int X, int Y,
			double scale, double sigma_scale, double quant, double ang_th,
			double log_eps, double density_th, int n_bins, int[][] reg_img,
			int[] reg_x, int[] reg_y) {
		image_double image;
		ntuple_list out = new ntuple_list(7);

		list_p = new coorlist();

		double[] return_value;
		image_double scaled_image, angles;
		image_char used;
		image_int region = null;
		// ArrayList<Point> list_p;
		// void * mem_p;
		rect rec = new rect();

		Point[] reg;
		// int reg_size = 0,
		int min_reg_size, i;
		int xsize, ysize;
		double rho, reg_angle = 0, prec, p, log_nfa, logNT;
		int ls_count = 0; /* line segments are numbered 1,2,3,... */

		/* check parameters检查参数的合理性 */
		if (img == null || X <= 0 || Y <= 0)
			error("invalid image input.");
		if (scale <= 0.0)
			error("'scale' value must be positive.");
		if (sigma_scale <= 0.0)
			error("'sigma_scale' value must be positive.");
		if (quant < 0.0)
			error("'quant' value must be positive.");
		if (ang_th <= 0.0 || ang_th >= 180.0)
			error("'ang_th' value must be in the range (0,180).");
		if (density_th < 0.0 || density_th > 1.0)
			error("'density_th' value must be in the range [0,1].");
		if (n_bins <= 0)
			error("'n_bins' value must be positive.");

		/* angle tolerance */
		prec = M_PI * ang_th / 180.0;
		p = ang_th / 180.0;
		rho = quant / Math.sin(prec); /* 梯度幅度阈值 */

		modgrad = new image_double(X, Y);

		/* 加载和缩放图像（如有必要）并计算每个像素的角度 */
		image = new image_double(X, Y, img);
		if (scale != 1.0) {
			scaled_image = gaussian_sampler(image, scale, sigma_scale);
			angles = ll_angle(scaled_image,/* list_p, */rho, /*
															 * &list_p, &mem_p,*
															 * &modgrad,
															 */(int) n_bins);
			// free_image_double(scaled_image);
		} else
			angles = ll_angle(image, /* list_p, */rho,/* list_p, &mem_p, &modgrad, */
					(int) n_bins);
		xsize = angles.xsize;
		ysize = angles.ysize;



       /*
         *测试次数 - NT
         *
         *理论上的测试次数是Np。（XY）^（5/2），其中X和Y是
         *图像的列数和行数。 Np对应于数字
         *考虑了角度精度。 正如过程'rect_improve'测试
         * 5倍减少角度精度，5次后
         *改进其他因素，有11个不同的精度值
         *可能经过测试。 因此，测试的数量是11 *（X * Y）^（5/2）
         *其对数值为log10（11）+ 5/2 *（log10（X）+ log10（Y））。
         */
		logNT = 5.0 * (Math.log10((double) xsize) + Math.log10((double) ysize))
				/ 2.0 + Math.log10(11.0);
		min_reg_size = (int) (-logNT / Math.log10(p)); /*
														 * minimal number of
														 * points in region that
														 * can give a meaningful
														 * event
														 */

		/* initialize some structures */
		// if( reg_img != NULL && reg_x != NULL && reg_y != NULL ) /* save
		// region data */
		region = new image_int(angles.xsize, angles.ysize, 0);
		used = new image_char(xsize, ysize, NOTUSED);
		reg = new Point[xsize * ysize];

		for (int z = 0; z < xsize * ysize; z++) {
			reg[z] = new Point();
		}
		// if( reg == NULL ) error("not enough memory!");

		//System.out.println("'ere1" + list_p + " " + list_p.next);
		/* search for line segments */
		for (; list_p != null; list_p = list_p.next) {
			// System.out.println("'ere1.5       "+(used.data[list_p.x +
			// list_p.y * used.xsize] == NOTUSED)+" "+(angles.data[list_p.x +
			// list_p.y * angles.xsize] != NOTDEF));
			if (used.data[list_p.x + list_p.y * used.xsize] == NOTUSED
					&& angles.data[list_p.x + list_p.y * angles.xsize] != NOTDEF)
			/*
			 * there is no risk of double comparison problems here because we
			 * are only interested in the exact NOTDEF value
			 */
			{
				// System.out.println("'ere2");

				/* find the region of connected point and ~equal angle */
				//System.out.println("attempting to grow " + list_p.x + " "
				//		+ list_p.y);
				region_grow(list_p.x, list_p.y, angles, reg, /*
															 * &reg_size,
															 * &reg_angle,
															 */used, prec);

				/* reject small regions */
				if (reg_size < min_reg_size) {
					//System.out.println("regsize " + reg_size + "   min: "
					//		+ min_reg_size);
					continue;
				}

				//System.out.println("LINE FOUND HERE");

				/* construct rectangular approximation for the region */
				region2rect(reg, reg_size, modgrad, /*reg_angle,*/ prec, p, rec);

				/*
				
               *检查矩形是否超出区域的最小密度
               *点。 如果不是，请尝试改善该地区。 矩形会
               *如果最后一个不符合最低标准，则被拒绝
               *密度条件。 这是对原始LSD的补充
               *算法发布在
               *“LSD：具有错误检测控制的快速线段检测器”
               * R. Grompone von Gioi，J. Jakubowicz，J.M. Morel和G.
               *兰德尔。 原始算法是用density_th =获得的
				 */
				if (!refine(reg, /*reg_size,*/ modgrad, /*reg_angle,*/ prec, p, rec,
						used, angles, density_th))
					continue;

				/* compute NFA value */
				log_nfa = rect_improve(rec, angles, logNT, log_eps);
				if (log_nfa <= log_eps)
					continue;

				/* A New Line Segment was found! */
				++ls_count; /* increase line segment counter */

				//System.out.println("LINE FOUND");
				/*
				*梯度是用2x2掩模计算的，它的值
                *对应于偏移量为（0.5,0.5）的点
                *应添加到输出。 坐标原点在
                *像素中心（0,0）。				 */
				rec.x1 += 0.5;
				rec.y1 += 0.5;
				rec.x2 += 0.5;
				rec.y2 += 0.5;

				/* scale the result values if a subsampling was performed */
				if (scale != 1.0) {
					rec.x1 /= scale;
					rec.y1 /= scale;
					rec.x2 /= scale;
					rec.y2 /= scale;
					rec.width /= scale;
				}

				/* add line segment found to output */
				
				//System.out.println(">>>>>>>> "+rec.x1+" "+ rec.y1+" "+ rec.x2+" "+ rec.y2);
				
				add_7tuple(out, rec.x1, rec.y1, rec.x2, rec.y2, rec.width,
						rec.p, log_nfa);

				/* add region number to 'region' image if needed */
				if (region != null)
					for (i = 0; i < reg_size; i++)
						region.data[reg[i].x + reg[i].y * region.xsize] = ls_count;
			}
		}

		 


	//	/ *返回结果* /
		// if（reg_img！= null && reg_x！= null && reg_y！= null）{
		/*
		* if（region == null）错误（“'region'应该是一个有效的图像。”）;
		* reg_img = region.data; 如果（region.xsize>（int）Integer.MAX_VALUE ||
		* region.xsize>（int）Integer.MAX_VALUE）
		*错误（“区域图像大到符合INT大小。”）; reg_x =（int）
		*（region.xsize）; reg_y =（int）（region.ysize）;
		*
		*
		 * 
		 * } if( out.size > ( int) Integer.MAX_VALUE )
		 * error("too many detections to fit in an INT.");n_out = (int)
		 * (out.size);
		 */		
		n_out = (int) (out.size);

		return_value = out.values;

		return return_value;
		// }

	}

	/*----------------------------------------------------------------------------*/
	/**
	 *LSD简单界面，带有比例和区域输出。
	 */
	
	double[] lsd_scale_region(double[] img, int X, int Y, double scale) {
		/* LSD 参数*/
		double sigma_scale = 0.3; /*
								 *用于高斯滤波器的西格玛计算为
                                 * sigma = sigma_scale / scale。
								 */
		double quant = 1.2; /*
						        *绑定到梯度上的量化误差
                                *规范。
							 */
		double ang_th = 37.5; /*以度为单位的梯度角公差。*/
		double log_eps = 0.0; /* 检测阈值：-log10（NFA）> log_eps */
		double density_th = 0.4; /*
								 * 区域的最小密度指向
								 *长方形。
								 */
		int n_bins = 1024; /*
							 *渐变伪排序中的箱数
                             *模量。
							 */

		return LineSegmentDetection(img, X, Y, scale, sigma_scale, quant,
				ang_th, log_eps, density_th, n_bins, null, null, null);

	}

	/*----------------------------------------------------------------------------*/
	/**
	 * LSD Simple Interface with Scale.//含有因子的调用
	 */
	double[] lsd_scale(double[] img, int X, int Y, double scale) {
		return lsd_scale_region(img, X, Y, scale/* ,null,null,null */);

	}

	public int n_out;

	/*----------------------------------------------------------------------------*/
	/**
	 * LSD Simple Interface.
	 */
	public double[] lsd(double[] img, int X, int Y) {
		/* LSD parameters */
		//double scale = 0.8; /* Scale the image by Gaussian filter to 'scale'. */
		double scale = 0.5;//固定因子
		return lsd_scale(img, X, Y, scale);//大小 ，高宽，因子
	}
	/*----------------------------------------------------------------------------*/

}
