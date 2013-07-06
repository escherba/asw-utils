/*
 * =====================================================================================
 *
 *       Filename:  AlignmentASW.hpp
 *
 *    Description:  Collection of routines for quality-aware alignment of Roche/454
 *                  reads. Class AlignmentASW implements asymmetric Smith-
 *                  Waterman-like algorithm with inverse scores (URL:
 *                  http://dx.doi.org/10.1101/gr.6468307)
 *
 *        Version:  1.0
 *        Created:  04/08/2011 15:55:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Eugene Scherba (es), escherba@bu.edu
 *        Company:  Boston University
 *
 * Copyright (c) 2011 Eugene Scherba
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * =====================================================================================
 */

#ifdef DEBUG
#define DEBUG
#endif

#include "boost/multi_array.hpp"
typedef uint32_t cigar_t;

template<typename T>
struct Matrix {
        typedef boost::multi_array<T, 2> type;
};


/*
 * =====================================================================================
 *        Class:  AlignmentASW
 *  Description:  
 * =====================================================================================
 */
class AlignmentASW
{
public:
        /* ====================  LIFECYCLE     ======================================= */
        AlignmentASW(int MATCH_PEN, int MISMATCH_PEN, int GAP_OPEN_EXTEND, int GAP_EXTEND);

        /* ====================  ACCESSORS     ======================================= */
        cigar_t* getCigar();
        size_t getCigarLength();
        int getAlignmentStart(int alstart_padded);
        int getScore();
        void checkScore(int score, uint32_t pos);
        void checkMatrix(cigar_t **other_matTra);
        void compareCigar(cigar_t* otherCig_p, cigar_t* otherCig_end);
        void print_cigar();
        void validateCigar();
#ifdef DEBUG
        void checkTrace();
        void checkTraceExtended();
        void validateCigarExtended();
        
#endif
/* ====================  MUTATORS      ======================================= */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  AlignmentASW
 *      Method:  AlignmentASW :: align
 * Description:  Fill out a traceback matrix and produce an alignment score
 *--------------------------------------------------------------------------------------
 */
        void align(const std::string &db,
		   const std::string &query,
		   const uint8_t *qual,
		   size_t clip_head,
		   size_t clip_tail);
/*
 *--------------------------------------------------------------------------------------
 *       Class:  AlignmentASW
 *      Method:  AlignmentASW :: trace
 * Description:  Create a CIGAR string from the traceback matrix
 *--------------------------------------------------------------------------------------
 */
        void trace();

        void softclip_trace();
/*
 *--------------------------------------------------------------------------------------
 *       Class:  AlignmentASW
 *      Method:  AlignmentASW :: append_softclip
 * Description:  extend CIGAR trace outer boundaries to previous clipping
 *--------------------------------------------------------------------------------------
 */
        void append_softclip();
/*
 *--------------------------------------------------------------------------------------
 *       Class:  AlignmentASW
 *      Method:  AlignmentASW :: append_hardclip
 * Description:  extend CIGAR trace outer boundaries to previous clipping
 *--------------------------------------------------------------------------------------
 */
        void append_hardclip(uint32_t clip_head, uint32_t clip_tail);

        void compact_trace();

        /* ====================  OPERATORS     ======================================= */

protected:
        /* ====================  DATA MEMBERS  ======================================= */

private:
        /* ====================  DATA MEMBERS  ======================================= */
        // Members below are assigned to during object initialization
        std::vector<int> match_penalty,
                         mismatch_penalty,
                         gopen_penalty,
                         gext_penalty;

        int GAP_OPEN_EXTEND,
            GAP_EXTEND;

        // Members below are assigned to in preparation for alignment
        const char *m_db,
	      	   *m_subdb,
                   *m_query,
		   *m_subquery;

        const uint8_t *m_qual,
	              *m_subqual;

	size_t m_db_len,
               m_subdb_len,
               m_query_len,
	       m_subquery_len;

        // Members below are assigned to during alignment
        std::vector<int> vecPen_m1_act,
                         vecPen_m_act,
                         vecIns_m1_act,
                         vecIns_m_act;
        std::vector<uint32_t> I_ext_m_act,
                              I_ext_m1_act;
        std::vector<int> *vecPen_m;

        std::vector<int>::const_iterator min_score_cell;
        Matrix<cigar_t>::type matTra;
#ifdef DEBUG
        Matrix<int>::type matPen,
                          matIns,
                          matDel;
#endif

        // Members below are assigned to during traceback
        int32_t offset;
        std::vector<cigar_t> rcigar;
        std::vector<cigar_t>::iterator fc5p;
        std::vector<cigar_t>::iterator fc3p;
}; /* -----  end of class AlignmentASW  ----- */
