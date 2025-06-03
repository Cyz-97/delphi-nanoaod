#ifndef SKELANA_ANALYSIS_HPP
#define SKELANA_ANALYSIS_HPP

#include <map>
#include <string>

#include "phdst_analysis.hpp"
#include "skelana.hpp"

namespace skelana
{
    class Analysis : public phdst::Analysis
    {
    private:
        typedef phdst::Analysis super;
    public:
        Analysis(const Analysis &) = delete;
        Analysis &operator=(const Analysis &) = delete;
        virtual ~Analysis();
        static Analysis *getInstance();
        void selectHadrons();
        void setOption(const std::string &option, int value);
        void setCut(const std::string &option, float value);
        void setIntCut(const std::string &option, int value);
        void printSkelanaSettings(std::ostream &os);
    protected:
        Analysis();
        virtual void user00();
        virtual int user01();
        virtual void user02();
        virtual void user99();
    private:
        bool checkHadron();
        bool checkRecordType();
        bool checkBadHPC();

        bool selectHadrons_;
        std::map<std::string, int> options_;
        inline static std::map<std::string, int*> optionMap_ =  {  
            {"IFLTRA", &pscflg_.ifltra},    
            {"IFLFIX", &pscflg_.iflfix},
            {"IFLRNQ", &pscflg_.iflrnq},
            {"IFLSTR", &pscflg_.iflstr},
            {"IFLJET", &pscflg_.ifljet},
            {"IFLSIM", &pscflg_.iflsim},
            {"IFLBSP", &pscflg_.iflbsp},
            {"IFLBTG", &pscflg_.iflbtg},
            {"IFLEMC", &pscflg_.iflemc},
            {"IFLHAC", &pscflg_.iflhac},
            {"IFLSTC", &pscflg_.iflstc},
            {"IFLELE", &pscflg_.iflele},
            {"IFLPHO", &pscflg_.iflpho},
            {"IFLMUO", &pscflg_.iflmuo},
            {"IFLHAD", &pscflg_.iflhad},
            {"IFLVDH", &pscflg_.iflvdh},
            {"IFLRV0", &pscflg_.iflrv0},
            {"IFLUTE", &pscflg_.iflute},
            {"IFLPHC", &pscflg_.iflphc},
            {"IFLVEC", &pscflg_.iflvec},
            {"IFLBHP", &pscflg_.iflbhp},
            {"IFLECL", &pscflg_.iflecl},
            {"IFLRVR", &pscflg_.iflrvr},
            {"IFLODR", &pscflg_.iflodr},
            {"IFLPVT", &pscflg_.iflpvt},
            {"IFLCUT", &pscflg_.iflcut},
            {"IFLVDR", &pscflg_.iflvdr},
            {"IFLFCT", &pscflg_.iflfct},
            {"IFLENR", &pscflg_.iflenr},

        };

        std::map<std::string, float> cuts_;
      inline static std::map<std::string, float*> cutMap_ =  {
        { "TRKMOM1", &pscutt_.trkmom[0] },
        { "TRKLEN1", &pscutt_.trklen[0] },
        { "TRKRPH1", &pscutt_.trkrph[0] },
        { "TRKZET1", &pscutt_.trkzet[0] },
        { "TRCCOS1", &pscutt_.trccos[0] },
        { "TRNCOS1", &pscutt_.trncos[0] },
        { "TRKMAX1", &pscutt_.trkmax[0] },
        { "TRKERR1", &pscutt_.trkerr[0] },
        { "RECCAL1", &pscutt_.reccal[0] },
        { "EHPC1", &pscutt_.ehpc[0] },
        { "EFEMC1", &pscutt_.efemc[0] },
        { "EHAC1", &pscutt_.ehac[0] },
        { "ESTIC1", &pscutt_.estic[0] },

            { "TRKMOM2", &pscutt_.trkmom[1] },
            { "TRKLEN2", &pscutt_.trklen[1] },
            { "TRKRPH2", &pscutt_.trkrph[1] },
            { "TRKZET2", &pscutt_.trkzet[1] },
            { "TRCCOS2", &pscutt_.trccos[1] },
            { "TRNCOS2", &pscutt_.trncos[1] },
            { "TRKMAX2", &pscutt_.trkmax[1] },
            { "TRKERR2", &pscutt_.trkerr[1] },
            { "RECCAL2", &pscutt_.reccal[1] },
            { "EHPC2", &pscutt_.ehpc[1] },
            { "EFEMC2", &pscutt_.efemc[1] },
            { "EHAC2", &pscutt_.ehac[1] },
            { "ESTIC2", &pscutt_.estic[1] },

            { "TRKMOM3", &pscutt_.trkmom[2] },
            { "TRKLEN3", &pscutt_.trklen[2] },
            { "TRKRPH3", &pscutt_.trkrph[2] },
            { "TRKZET3", &pscutt_.trkzet[2] },
            { "TRCCOS3", &pscutt_.trccos[2] },
            { "TRNCOS3", &pscutt_.trncos[2] },
            { "TRKMAX3", &pscutt_.trkmax[2] },
            { "TRKERR3", &pscutt_.trkerr[2] },
            { "RECCAL3", &pscutt_.reccal[2] },
            { "EHPC3", &pscutt_.ehpc[2] },
            { "EFEMC3", &pscutt_.efemc[2] },
            { "EHAC3", &pscutt_.ehac[2] },
            { "ESTIC3", &pscutt_.estic[2] },
        };

        std::map<std::string, int> intCuts_;
      inline static std::map<std::string, int*> intCutMap_ =  {
            { "VDONLY1", &pscutt_.vdonly[0] },
            { "IDVDWZ1", &pscutt_.idvdwz[0] },
            { "IHADRJ1", &pscutt_.ihadrj[0] },
            { "ISTOEL1", &pscutt_.istoel[0] },

            { "VDONLY2", &pscutt_.vdonly[1] },
            { "IDVDWZ2", &pscutt_.idvdwz[1] },
            { "IHADRJ2", &pscutt_.ihadrj[1] },
            { "ISTOEL2", &pscutt_.istoel[1] },

            { "VDONLY3", &pscutt_.vdonly[2] },
            { "IDVDWZ3", &pscutt_.idvdwz[2] },
            { "IHADRJ3", &pscutt_.ihadrj[2] },
            { "ISTOEL3", &pscutt_.istoel[2] }
      };
    };
};

#endif // SKELANA_ANALYSIS_HPP
