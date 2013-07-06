/*
 * =====================================================================================
 *
 *       Filename:  AlignmentASW.cpp
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
 * =====================================================================================
 */

// Standard C++ includes
#include <iostream>
//#include <map>
#include <cassert>
#include <algorithm>
#include <iomanip>
//#include <fstream>
//#include <iterator>
#include <vector>
//#include "bam.h"

#include "AlignmentASW.hpp"

/*
#ifdef NDEBUG
#define assert(x) \
{\
        if (!(x)) {\
                abort();\
        }\
}
#endif
*/

/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

/*
  CIGAR operations.
 */
/*! @abstract CIGAR: match */
#define BAM_CMATCH      0
/*! @abstract CIGAR: insertion in the read/donor, deletion in reference */
#define BAM_CINS        1
/*! @abstract CIGAR: deletion in the read/donor, insertion in reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: clip on the read with clipped sequence present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: padding */
#define BAM_CPAD        6
#define BAM_CSEQ_MATCH        7
#define BAM_CSEQ_MISMATCH 8

// Sanger PHRED scores range from 0 to 93
#define MAX_PHRED 94

#define AMBIGUOUS_BASE 'N'

/*
 * CIGAR operation codes (http://samtools.sourceforge.net/SAM-1.3.pdf):
 *
 * M - match or mismatch
 * I - insertion in the read/donor, deletion in reference
 * D - deletion in the read/donor, insertion in reference
 * N - skipped region from reference
 * S - soft-clip in the read
 * H - hard-clip in the read
 * P - padding (silent deletion from padded reference)
 * = - match
 * X - mismatch
 *                           0    1    2    3    4    5    6    7    8 */
const char cigar_chars[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};



template<typename T>
struct Submatrix {
        //typedef typename Matrix<T>::type::subarray<1>::type type;
        typedef typename Matrix<T>::type::template subarray<1>::type type;
};

template<typename T>
struct ConstSubmatrix {
        //typedef typename Matrix<T>::type::subarray<1>::type type;
        typedef typename Matrix<T>::type::template const_subarray<1>::type type;
};


template <typename T, typename V>
inline void print_vector(std::ostream& f, const T &vec, int width=0)
{
        typename T::const_iterator cii;
        for (cii = vec.begin(); cii != vec.end(); ++cii) {
                f << std::setw(width) << static_cast<V>(*cii);
        }
        f << std::endl;
}

template <typename T>
void print_align_matrix(std::ostream& f,
                        const std::string &db, 
                        const std::string &query, 
                        const uint8_t *qual, 
                        const typename Matrix<T>::type &A, 
                        int width=0)
{
        typedef typename ConstSubmatrix<T>::type ConstSubmatrixT;
        typename Matrix<T>::type::const_iterator i;
        typename ConstSubmatrixT::const_iterator j;
        std::string::const_iterator vii;
        const uint8_t *uii;

        i = A.begin();
        f << std::setw(width) << "" 
          << std::setw(width) << "" 
          << std::setw(width) << "";
        unsigned int col_ind = 0u;
        for (j = i->begin(); j != i->end(); ++j, ++col_ind) {
                f << std::setw(width) << col_ind;
        }
        f << std::endl;

        f << std::setw(width) << "" 
          << std::setw(width) << "" 
          << std::setw(width) << "" 
          << std::setw(width) << "-";
        print_vector<std::string, char>(f, db, width);

        i = A.begin();
        f << std::setw(width) << 0u 
          << std::setw(width) << "-" 
          << std::setw(width) << "-";
        print_vector<ConstSubmatrixT, int>(f, *i, width);
        ++i;

        // The rest of the matrix
        unsigned int k = 1u;
        for (vii = query.begin(), uii = qual;
             i != A.end();
             ++i, ++vii, ++uii, ++k)
        {
                f << std::setw(width) << k 
                  << std::setw(width) << static_cast<int>(*uii) 
                  << std::setw(width) << *vii;
                print_vector<ConstSubmatrixT, int>(f, *i, width);
        }
        f << std::endl;
}

inline bool match_ref(char cRef, char cQuery) {
        // Asymmetric match: avoid regions containing Ns on the reference
        // but do not penalize for Ns on the query
        return (cRef == cQuery || cQuery == AMBIGUOUS_BASE);
}

AlignmentASW::AlignmentASW(int MATCH_PEN, int MISMATCH_PEN, int pGAP_OPEN_EXTEND, int pGAP_EXTEND) :
        match_penalty(std::vector<int>(MAX_PHRED)),
        mismatch_penalty(std::vector<int>(MAX_PHRED)),
        gopen_penalty(std::vector<int>(MAX_PHRED)),
        gext_penalty(std::vector<int>(MAX_PHRED)),
        GAP_OPEN_EXTEND(pGAP_OPEN_EXTEND),
        GAP_EXTEND(pGAP_EXTEND),
        offset(0)
{
        // Initialize penalty look-up vectors
        const double qN = -10.0 * log10(0.75); // P(error | N) = 0.75
        for (std::vector<int>::size_type i = 0u; i < MAX_PHRED; ++i) {
                double weight = 1.0 - pow(10.0, -(static_cast<double>(i) + qN)/10.0);
                match_penalty[i] = 10 + static_cast<int>(round(weight * static_cast<double>(MATCH_PEN)));
                mismatch_penalty[i] = 10 + static_cast<int>(round(weight * static_cast<double>(MISMATCH_PEN)));
                gopen_penalty[i] = 10 + static_cast<int>(round(weight * static_cast<double>(GAP_OPEN_EXTEND)));
                gext_penalty[i] = 10 + static_cast<int>(round(weight * static_cast<double>(GAP_EXTEND)));
        }
        //print_vector<std::vector<int>, int>(std::cout, match_penalty, 5);
        //print_vector<std::vector<int>, int>(std::cout, mismatch_penalty, 5);
        //print_vector<std::vector<int>, int>(std::cout, gopen_penalty, 5);
        //print_vector<std::vector<int>, int>(std::cout, gext_penalty, 5);
}

/* ====================  ACCESSORS     ======================================= */
cigar_t* AlignmentASW::getCigar() {
        return &(*fc5p);
}
size_t AlignmentASW::getCigarLength() {
#ifdef DEBUG
        assert(fc3p > fc5p);
#endif
        return fc3p - fc5p;
}
int AlignmentASW::getAlignmentStart(int alstart_padded) {
        assert(m_subdb >= m_db);
        return std::max<int>(0, alstart_padded) + offset + (m_subdb - m_db);
}
int AlignmentASW::getScore() {
        return *min_score_cell;
}
void AlignmentASW::checkScore(int score, uint32_t pos) {
        assert (score == *min_score_cell);
        assert((min_score_cell - (*vecPen_m).begin()) == static_cast<int>(pos));
}
void AlignmentASW::checkMatrix(cigar_t **other_matTra) {
        Matrix<cigar_t>::type::index i, j;
        Matrix<cigar_t>::type::index imax = matTra.size();
        Submatrix<cigar_t>::type rowTra = matTra[0];
        Matrix<cigar_t>::type::index jmax = rowTra.size();
        for (i = 0; i < imax; ++i) {
                for (j = 0; j < jmax; ++j) {
                        assert(matTra[i][j] == other_matTra[i][j]);
                }
        }
}
void AlignmentASW::compareCigar(cigar_t* otherCig_p, cigar_t* otherCig_end) {
        assert(fc3p - fc5p == otherCig_end - otherCig_p);
        std::vector<cigar_t>::iterator c = fc5p;
        for (; c != fc3p; ++c, ++otherCig_p) {
                assert(*c == *otherCig_p);
        }
}
void AlignmentASW::print_cigar() {
        std::vector<cigar_t>::iterator c;
        for (c = fc5p; c != fc3p; ++c) {
                cigar_t cigar = *c;
                std::cout << (cigar >> BAM_CIGAR_SHIFT) 
                          << cigar_chars[cigar & BAM_CIGAR_MASK] 
                          << ' ';
        }
        std::cout << std::endl;
}
void AlignmentASW::validateCigar() {
        std::vector<cigar_t>::iterator c;
        size_t exp_len = 0u;
        if(fc3p - fc5p <= 0) {
                abort();
        }
        for (c = fc5p; c != fc3p; ++c) {
                cigar_t cigar = *c;
                int state = cigar & BAM_CIGAR_MASK;
                uint32_t ops = cigar >> BAM_CIGAR_SHIFT;
                switch(state) {
                case BAM_CMATCH:
                case BAM_CINS:
                case BAM_CSOFT_CLIP:
                        exp_len += ops;
                        break;
                case BAM_CDEL:
                        break;
                default:
                        std::cerr << "Invalid CIGAR character " << state 
                                  << " at position " << c - fc5p 
                                  << std::endl;
                        print_cigar();
                        abort();
                        break;
                }
        }
        assert(exp_len == matTra.size() - 1u);
}
#ifdef DEBUG
void AlignmentASW::checkTrace() {
        std::vector<cigar_t>::iterator c = fc5p;
        const char *db = m_subdb + offset,
                   *query = m_subquery;
        const uint8_t *qual = m_subqual;

        Matrix<int>::type::index m1 = 0,
                                 n1 = offset;
        int score = 0;
        for (; c != rcigar.end(); ++c) { // Loop over the CIGAR string
                int l = *c >> BAM_CIGAR_SHIFT;
                assert(l > 0);
                int op = *c & BAM_CIGAR_MASK;
                assert(op >= 0);

                switch (op) {
                case BAM_CMATCH:
                        for (int i = 0; i < l; ++i, ++db, ++query, ++qual, ++m1, ++n1) {
                                score += (match_ref(*db, *query) 
                                     ? match_penalty[static_cast<int>(*qual)] 
                                     : mismatch_penalty[static_cast<int>(*qual)]);
                        }
                        break;
                case BAM_CDEL:
                        score += GAP_OPEN_EXTEND;
                        ++db, ++n1;
                        for (int i = 1; i < l; ++i, ++db, ++n1) {
                                score += GAP_EXTEND;
                        }
                        break;
                case BAM_CINS:
                        score += gopen_penalty[static_cast<int>(*qual)];
                        ++query, ++qual, ++m1;
                        for (int i = 1; i < l; ++i, ++query, ++qual, ++m1) {
                                score += gext_penalty[static_cast<int>(*qual)];
                        }
                        break;
                default:
                        std::cout << "ERROR: Unknown CIGAR operation " << op 
                                  << std::endl;
                        abort();
                        break;
                }
        }
        assert(*min_score_cell == score);
}

void AlignmentASW::checkTraceExtended() {
        std::vector<cigar_t>::iterator c;
        const char *db = m_subdb + offset,
                   *query = m_subquery;
        const uint8_t *qual = m_subqual;

        Matrix<int>::type::index m1 = 0,
                                 n1 = offset;
        int score = 0;
        for (c = fc5p; c != fc3p; ++c) { // Loop over the CIGAR string
                int l = *c >> BAM_CIGAR_SHIFT;
                assert(l > 0);
                int op = *c & BAM_CIGAR_MASK;
                assert(op >= 0);

                switch (op) {
                case BAM_CSEQ_MATCH:
                        for (int i = 0; i < l; ++i) {
                                score += match_penalty[static_cast<int>(*qual)];
                                ++db, ++query, ++qual, ++m1, ++n1;
                        }
                        break;
                case BAM_CSEQ_MISMATCH:
                        for (int i = 0; i < l; ++i) {
                                score += mismatch_penalty[static_cast<int>(*qual)];
                                ++db, ++query, ++qual, ++m1, ++n1;
                        }
                        break;
                case BAM_CDEL:
                        score += GAP_OPEN_EXTEND;
                        ++db, ++n1;
                        for (int i = 1; i < l; ++i) {
                                score += GAP_EXTEND;
                                ++db, ++n1;
                        }
                        break;
                case BAM_CINS:
                        score += gopen_penalty[static_cast<int>(*qual)];
                        ++query, ++qual, ++m1;
                        for (int i = 1; i < l; ++i) {
                                score += gext_penalty[static_cast<int>(*qual)];
                                ++query, ++qual, ++m1;
                        }
                        break;
                default:
                        std::cout << "ERROR: Unknown CIGAR operation " << op 
                                  << std::endl;
                        abort();
                        break;
                }
        }
        assert(*min_score_cell == score);
}

void AlignmentASW::validateCigarExtended() {
        std::vector<cigar_t>::iterator c;
        size_t exp_len = 0u;
        assert(fc3p - fc5p > 0);
        for (c = fc5p; c != fc3p; ++c) {
                cigar_t cigar = *c;
                int state = cigar & BAM_CIGAR_MASK;
                uint32_t ops = cigar >> BAM_CIGAR_SHIFT;
                switch(state) {
                case BAM_CSEQ_MATCH:
                case BAM_CSEQ_MISMATCH:
                case BAM_CINS:
                        exp_len += ops;
                        break;
                case BAM_CDEL:
                        break;
                default:
                        std::cerr << "Invalid CIGAR character " << state 
                                  << " at position " << c - fc5p 
                                  << std::endl;
                        print_cigar();
                        abort();
                        break;
                }
        }
        assert(exp_len = matTra.size() - 1);
}

#endif
/* ====================  MUTATORS      ======================================= */

/*
*--------------------------------------------------------------------------------------
*       Class:  AlignmentASW
*      Method:  AlignmentASW :: align
* Description:  Fill out a traceback matrix and produce an alignment score
*--------------------------------------------------------------------------------------
*/
void AlignmentASW::align(const std::string &db,
           const std::string &query,
           const uint8_t *qual,
           size_t clip_head,
           size_t clip_tail)
{
        m_db_len = db.size(),
        m_query_len = query.size();
        m_subdb_len = m_db_len - clip_head - clip_tail,
        m_subquery_len = m_query_len - clip_head - clip_tail;

        Matrix<int>::type::index m, m1, n, n1;

        m_db = db.c_str();
        m_query = query.c_str();
        m_qual = qual;
        
        m_subdb = m_db + clip_head;
        m_subquery = m_query + clip_head;
        m_subqual = m_qual + clip_head;


        // vectors below store scores
        vecPen_m1_act.resize(m_subdb_len + 1);
        vecPen_m_act.resize(m_subdb_len + 1);
        vecIns_m1_act.resize(m_subdb_len + 1);
        vecIns_m_act.resize(m_subdb_len + 1);

        // vectors below store extension counts
        I_ext_m_act.resize(m_subdb_len + 1);
        I_ext_m1_act.resize(m_subdb_len + 1);

        // actual edit matrix
        matTra.resize(boost::extents[m_subquery_len + 1][m_subdb_len + 1]);
#ifdef DEBUG
        matPen.resize(boost::extents[m_subquery_len + 1][m_subdb_len + 1]);
        matIns.resize(boost::extents[m_subquery_len + 1][m_subdb_len + 1]);
        matDel.resize(boost::extents[m_subquery_len + 1][m_subdb_len + 1]);
#endif

        std::vector<int> *vecPen_m1 = &vecPen_m1_act,
                         *vecIns_m1 = &vecIns_m1_act,
                         *vecIns_m = &vecIns_m_act;
        vecPen_m = &vecPen_m_act;

        std::vector<uint32_t> *I_ext_m = &I_ext_m_act,
                              *I_ext_m1 = &I_ext_m1_act;

        // Initialize first row
        //
        // qq - quality in the query at position m
        unsigned int qq = static_cast<unsigned int>(m_subqual[0]);
        int gopen_true_pen = gopen_penalty[qq] - gext_penalty[qq];

        vecPen_m_act[0] = 0;
        vecIns_m_act[0] = vecPen_m_act[0] + gopen_true_pen;
        I_ext_m_act[0] = 0;
        int storedDel_score = vecPen_m_act[0] + (GAP_OPEN_EXTEND - GAP_EXTEND);
#ifdef DEBUG
        matPen[0][0] = vecPen_m_act[0];
        matDel[0][0] = storedDel_score;
        matIns[0][0] = vecIns_m_act[0];
#endif
        Matrix<cigar_t>::type::iterator rowTra = matTra.begin();
        Submatrix<cigar_t>::type::iterator cellTra = (*rowTra).begin();

        *cellTra = (0 << BAM_CIGAR_SHIFT) | BAM_CSEQ_MATCH;
        ++cellTra;
        for (n = 0, n1 = 1; n < static_cast<int>(m_subdb_len); ++n, ++n1) {
                vecPen_m_act[n1] = 0;
                vecIns_m_act[n1] = vecPen_m_act[n1] + gopen_true_pen;
                //                       ^
                //                       zero for semi-global alignment
                I_ext_m_act[n1] = 0;

                *cellTra = (n1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
#ifdef DEBUG
                matPen[0][n1] = vecPen_m_act[n1];
                matDel[0][n1] = storedDel_score;
                matIns[0][n1] = vecIns_m_act[n1];
#endif
                ++cellTra;
        }
        ++rowTra;

        for (m = 0, m1 = 1; m < static_cast<int>(m_subquery_len); ++m, ++m1) {
                // cq - character in the query at position m
                char cq = m_subquery[m];
                // qq - quality in the query at position m
                unsigned int qq = static_cast<unsigned int>(m_subqual[m]);
                int match_pen = match_penalty[qq],
                    mismatch_pen = mismatch_penalty[qq],
                    gopen_pen = gopen_penalty[qq],
                    gext_pen = gext_penalty[qq];

                cellTra = (*rowTra).begin();
                std::vector<int>::iterator iPen_m = (*vecPen_m).begin();
                std::vector<int>::iterator iPen_m1 = (*vecPen_m1).begin();
                std::vector<int>::iterator iIns_m = (*vecIns_m).begin();
                std::vector<int>::iterator iIns_m1 = (*vecIns_m1).begin();
                std::vector<uint32_t>::iterator iExt_m = (*I_ext_m).begin();
                std::vector<uint32_t>::iterator iExt_m1 = (*I_ext_m1).begin();

                uint32_t cI;
                int wI_extend;
                *iIns_m1 = wI_extend = *iIns_m + gext_pen;
                *iExt_m1 = cI = *iExt_m + 1u;
                *cellTra = (cI << BAM_CIGAR_SHIFT) | BAM_CINS;

                int Pen_m1_n;
                *iPen_m1 = Pen_m1_n = wI_extend;
                storedDel_score = wI_extend + (GAP_OPEN_EXTEND - GAP_EXTEND);
                //                  ^
                //                  equal to *iPen_m1
#ifdef DEBUG
                Submatrix<int>::type rowPen = matPen[m1];
                Submatrix<int>::type rowDel = matDel[m1];
                Submatrix<int>::type rowIns = matIns[m1];
                rowPen[0] = *iPen_m1;
                rowIns[0] = *iIns_m1;
                rowDel[0] = storedDel_score;
#endif
                int Pen_m_n = *iPen_m;

                ++cellTra;
                ++iPen_m, ++iPen_m1;
                ++iIns_m, ++iIns_m1;
                ++iExt_m, ++iExt_m1;

                uint32_t cD = 0u;
                for (n = 0, n1 = 1; n < static_cast<int>(m_subdb_len); ++n, ++n1) {
                        int wD, wI, wM;
                        uint32_t cI;
                        bool is_seq_match = match_ref(m_subdb[n], cq);

                        int wD_open = Pen_m1_n + GAP_OPEN_EXTEND;
                        int wD_extend = storedDel_score + GAP_EXTEND;

                        int wI_open = *iPen_m + gopen_pen;
                        int wI_extend = *iIns_m + gext_pen;

                        // given equal scores, prefer extending 
                        // existing gaps to opening new ones
                        if (wD_open < wD_extend) {
                                storedDel_score = wD = wD_open;
                                cD = 1u;
                        } else {
                                storedDel_score = wD = wD_extend;
                                ++cD;
                        }
                        if (wI_open < wI_extend) {
                                *iIns_m1 = wI = wI_open;
                                *iExt_m1 = cI = 1u;
                        } else {
                                *iIns_m1 = wI = wI_extend;
                                *iExt_m1 = cI = *iExt_m + 1u;
                        }

                        int mstate;
                        if (is_seq_match) {
                                wM = Pen_m_n + match_pen;
                                mstate = BAM_CSEQ_MATCH;
                        } else {
                                wM = Pen_m_n + mismatch_pen;
                                mstate = BAM_CSEQ_MISMATCH;
                        }

                        // Order of preference: M, I, D
                        if (wI < wM) {
                                // either insertion or deletion
                                if (wD < wI) {
                                        // deletion
                                        *cellTra = (cD << BAM_CIGAR_SHIFT) | BAM_CDEL;
                                        *iPen_m1 = wD;
                                } else {
                                        // insertion
                                        *cellTra = (cI << BAM_CIGAR_SHIFT) | BAM_CINS;
                                        *iPen_m1 = wI;
                                }
                        } else if (wD < wM) {
                                // deletion
                                *cellTra = (cD << BAM_CIGAR_SHIFT) | BAM_CDEL;
                                *iPen_m1 = wD;
                        } else {
                                // either match or mismatch
                                *cellTra = (1u << BAM_CIGAR_SHIFT) | mstate;
                                *iPen_m1 = wM;
                        }
#ifdef DEBUG
                        rowDel[n1] = wD;
                        rowIns[n1] = wI;
                        //rowPen[n1] = (*vecPen_m1)[n1];
                        rowPen[n1] = *iPen_m1;
#endif
                        Pen_m1_n = *iPen_m1;
                        Pen_m_n = *iPen_m;

                        ++cellTra;
                        ++iPen_m, ++iPen_m1;
                        ++iIns_m, ++iIns_m1;
                        ++iExt_m, ++iExt_m1;
                }
                std::swap(vecIns_m1, vecIns_m);
                std::swap(vecPen_m1, vecPen_m);

                std::swap(I_ext_m1, I_ext_m);
                // At this point, vecIns_m and vecPen_m point to their
                // corresponding last rows, and vecIns_m1 and vecPen_m1 
                // point to penultimate rows
                ++rowTra;
        }

        // Find minimum penalty in final row
        min_score_cell = std::min_element((*vecPen_m).begin(), (*vecPen_m).end());
#ifdef DEBUG
        //print_align_matrix<int>(std::cout, m_subdb, m_subquery, m_subqual, matPen, 5);
        //print_align_matrix<int>(std::cout, m_subdb, m_subquery, m_subqual, matIns, 5);
        //print_align_matrix<int>(std::cout, m_subdb, m_subquery, m_subqual, matDel, 5);
        //print_align_matrix<cigar_t>(std::cout, m_subdb, m_subquery, m_subqual, matTra, 5);
#endif
}

/*
*--------------------------------------------------------------------------------------
*       Class:  AlignmentASW
*      Method:  AlignmentASW :: trace
* Description:  Create a CIGAR string from the traceback matrix
*--------------------------------------------------------------------------------------
*/
void AlignmentASW::trace() {
        assert(matTra.size() == m_subquery_len + 1);
        rcigar.resize(m_subquery_len + 2u);
        Matrix<int>::type::index m1 = m_subquery_len,
                                 n1 = min_score_cell - (*vecPen_m).begin();

        // fill out cigar string
        cigar_t cigar = matTra[m1][n1];
        // z - length of CIGAR operation
        uint32_t z = cigar >> BAM_CIGAR_SHIFT;
        // state - CIGAR operation
        uint32_t state = cigar & BAM_CIGAR_MASK;

        std::vector<cigar_t>::reverse_iterator rc = rcigar.rbegin();
        ++rc; // skip last element

        std::vector<cigar_t>::reverse_iterator start_riter = rc;

        uint32_t num_matches = 0u;

        while (m1 > 0) {
        switch (state) {
        case BAM_CSEQ_MATCH:
                num_matches = 0;
                do {
                        // simply accumulate num_matches
                        num_matches += z;
                        m1 -= z, n1 -= z;
                        cigar = matTra[m1][n1];
                        z = cigar >> BAM_CIGAR_SHIFT;
                        state = cigar & BAM_CIGAR_MASK;
                } while (state == BAM_CSEQ_MATCH && m1 > 0);
                *rc = (num_matches << BAM_CIGAR_SHIFT) | BAM_CSEQ_MATCH;
                ++rc;
                break;
        case BAM_CSEQ_MISMATCH:
                num_matches = 0;
                do {
                        // simply accumulate num_matches
                        num_matches += z;
                        m1 -= z, n1 -= z;
                        cigar = matTra[m1][n1];
                        z = cigar >> BAM_CIGAR_SHIFT;
                        state = cigar & BAM_CIGAR_MASK;
                } while (state == BAM_CSEQ_MISMATCH && m1 > 0);
                *rc = (num_matches << BAM_CIGAR_SHIFT) | BAM_CSEQ_MISMATCH;
                ++rc;
                break;
        case BAM_CDEL:
                *rc = cigar;
                ++rc;
                n1 -= z;
                cigar = matTra[m1][n1];
                z = cigar >> BAM_CIGAR_SHIFT;
                state = cigar & BAM_CIGAR_MASK;
                break;
        case BAM_CINS:
                *rc = cigar;
                ++rc;
                m1 -= z;
                cigar = matTra[m1][n1];
                z = cigar >> BAM_CIGAR_SHIFT;
                state = cigar & BAM_CIGAR_MASK;
                break;
        default:
                std::cout << "ERROR (trace): unknown CIGAR operation " << state 
                          << std::endl;
                abort();
                break;
        }}
        offset = n1;
        fc5p = rc.base();
        fc3p = rcigar.end();
        --fc3p;

#ifdef DEBUG
        //print_cigar();
        validateCigarExtended();
        checkTraceExtended();
#endif
}

void AlignmentASW::softclip_trace()
{
        // backtrack to the last match: shift fc3p
        unsigned int soft_clip_3p = 0u;
        std::vector<cigar_t>::reverse_iterator rc, rc3p;
        rc3p = std::vector<cigar_t>::reverse_iterator(fc3p); //rcigar.rbegin()
        rc = std::vector<cigar_t>::reverse_iterator(fc5p);

        for (; rc3p != rc; ++rc3p) {
                cigar_t cigar = *rc3p;
                uint32_t state = cigar & BAM_CIGAR_MASK;
                if (state == BAM_CSEQ_MATCH) {
                        break;
                } else if (state != BAM_CDEL) {
                        // BAM_CSEQ_MISMATCH or BAM_CINS
                        soft_clip_3p += (cigar >> BAM_CIGAR_SHIFT);
                }
        }
        if (soft_clip_3p > 0u) {
                --rc3p;
                *rc3p = (soft_clip_3p << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
        }
        fc3p = rc3p.base();

        // track to the first match: shift fc5p
        unsigned int soft_clip_5p = 0u;
        for (; fc5p != fc3p; ++fc5p) {
                cigar_t cigar = *fc5p;
                uint32_t state = cigar & BAM_CIGAR_MASK;
                if (state == BAM_CSEQ_MATCH) {
                        break;
                } else {
                        int op_len = cigar >> BAM_CIGAR_SHIFT;
                        if (state != BAM_CDEL) {
                                // BAM_CSEQ_MISMATCH or BAM_CINS
                                soft_clip_5p += op_len;
                        }
                        if (state != BAM_CINS) {
                                // BAM_CDEL, BAM_CSEQ_MISMATCH
                                offset += op_len;
                        }
                }
        }
        if (soft_clip_5p > 0u) {
                --fc5p;
                *fc5p = (soft_clip_5p << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP;
        }
}

/*
*--------------------------------------------------------------------------------------
*       Class:  AlignmentASW
*      Method:  AlignmentASW :: append_softclip
* Description:  extend CIGAR trace outer boundaries to previous clipping
*--------------------------------------------------------------------------------------
*/
void AlignmentASW::append_softclip()
{
        assert(m_subquery >= m_query);
        assert(m_subdb >= m_db);
        uint32_t clip_head = m_subquery - m_query;
        if (clip_head > 0u) {
                /* clipped beginnning */
                cigar_t cigar = *fc5p;
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CSOFT_CLIP) {
                        /* extend clipping */
                        *fc5p = ((clip_head + z) << BAM_CIGAR_SHIFT)
                                       | BAM_CSOFT_CLIP;
                } else if (state == BAM_CSEQ_MATCH || state == BAM_CMATCH) {
                        /* try to contract clipping */
                        uint32_t match_add = 0u;
                        const char *subquery = m_subquery,
                                   *subdb = m_subdb + offset;
                        while (clip_head > 0u && *--subquery == *--subdb) {
                                ++match_add;
                                --clip_head;
                        }
                        if (match_add > 0u) {
                                *fc5p = ((z + match_add) << BAM_CIGAR_SHIFT)
                                               | state;
                                offset -= match_add;
                        }
                        if (clip_head > 0u) {
                                *--fc5p = (clip_head << BAM_CIGAR_SHIFT)
                                                   | BAM_CSOFT_CLIP;
                        }
                } else {
                        /* add clipping */
                        *--fc5p = (clip_head << BAM_CIGAR_SHIFT)
                                               | BAM_CSOFT_CLIP;
                }
        }
        uint32_t clip_tail = (m_query + m_query_len) - (m_subquery + m_subquery_len);
        if (clip_tail > 0u) {
                /* clipped end */
                cigar_t cigar = *(fc3p - 1u);
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CSOFT_CLIP) {
                        /* extend clipping */
                        *(fc3p - 1u) = ((clip_tail + z) << BAM_CIGAR_SHIFT)
                                    | BAM_CSOFT_CLIP;
                } else if (state == BAM_CSEQ_MATCH || state == BAM_CMATCH) {
                        /* try to contract clipping */
                        uint32_t match_add = 0u;
                        const char *subquery = m_subquery + m_subquery_len - 1u,
                                   *subdb = m_subdb + offset + m_subdb_len - 1u;
                        while (clip_tail > 0u && *++subquery == *++subdb) {
                                ++match_add;
                                --clip_tail;
                        }
                        if (match_add > 0u) {
                                *(fc3p - 1u)
                                        = ((z + match_add) << BAM_CIGAR_SHIFT) | state;
                        }
                        if (clip_tail > 0u) {
                                *fc3p = (clip_tail << BAM_CIGAR_SHIFT)
                                                   | BAM_CSOFT_CLIP;
                                ++fc3p;
                        }
                } else {
                        /* add clipping */
                        *fc3p = (clip_tail << BAM_CIGAR_SHIFT)
                                               | BAM_CSOFT_CLIP;
                        ++fc3p;
                }
        }
}
/*
*--------------------------------------------------------------------------------------
*       Class:  AlignmentASW
*      Method:  AlignmentASW :: append_hardclip
* Description:  extend CIGAR trace outer boundaries to previous clipping
*--------------------------------------------------------------------------------------
*/
void AlignmentASW::append_hardclip(uint32_t clip_head, uint32_t clip_tail)
{
        if (clip_head > 0u) {
                /* clipped beginnning */
                cigar_t cigar = *fc5p;
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CHARD_CLIP) {
                        /* extend clipping */
                        *fc5p = ((clip_head + z) << BAM_CIGAR_SHIFT)
                                       | BAM_CHARD_CLIP;
                } else {
                        /* add clipping */
                        *--fc5p = (clip_head << BAM_CIGAR_SHIFT)
                                               | BAM_CHARD_CLIP;
                }
        }
        if (clip_tail > 0u) {
                /* clipped end */
                cigar_t cigar = *(fc3p - 1u);
                uint32_t state = cigar & BAM_CIGAR_MASK,
                         z = cigar >> BAM_CIGAR_SHIFT;
                if (state == BAM_CHARD_CLIP) {
                        /* extend clipping */
                        *(fc3p - 1u) = ((clip_tail + z) << BAM_CIGAR_SHIFT)
                                    | BAM_CHARD_CLIP;
                } else {
                        /* add clipping */
                        *fc3p = (clip_tail << BAM_CIGAR_SHIFT)
                                               | BAM_CHARD_CLIP;
                        ++fc3p;
                }
        }
}

void AlignmentASW::compact_trace()
{
        // collapse the CIGAR string by treating matches 
        // and mismatches as the same state
        std::vector<cigar_t>::reverse_iterator rbucket, start_riter, rc;
        rbucket = start_riter 
                = std::vector<cigar_t>::reverse_iterator(fc3p); //rcigar.rbegin()
        rc = std::vector<cigar_t>::reverse_iterator(fc5p);

        while (rbucket != rc)
        {
                cigar_t cigar;
                uint32_t num_matches = 0u;
                while (true) {
                        cigar = *rbucket;
                        uint32_t state = cigar & BAM_CIGAR_MASK;
                        ++rbucket;
                        if (state == BAM_CSEQ_MATCH 
                         || state == BAM_CSEQ_MISMATCH) {
                                num_matches += (cigar >> BAM_CIGAR_SHIFT);
                                if (rbucket == rc) {
                                        cigar = (num_matches << BAM_CIGAR_SHIFT)
                                              | BAM_CMATCH;
                                        break;
                                }
                        } else if (num_matches > 0u) {
                                *start_riter = (num_matches << BAM_CIGAR_SHIFT)
                                             | BAM_CMATCH;
                                ++start_riter;
                                break;
                        } else {
                                break;
                        }
                }
                *start_riter = cigar;
                ++start_riter;
        }
        fc5p = start_riter.base();
#ifdef NDEBUG
        // release mode
        //validateCigar();
#else
        // debug mode
        //assert(fc3p >= fc5p);
        //print_cigar();
        //checkTrace();
#endif
}

