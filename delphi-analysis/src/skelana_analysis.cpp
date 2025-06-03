#include <iostream>
#include <iomanip>
#include <set>
#include <string>


#include "skelana_analysis.hpp"
#include "skelana.hpp"
#include "phdst.hpp"

namespace skelana
{
    Analysis::Analysis() : selectHadrons_(false) {};

    Analysis::~Analysis() {};

    Analysis *Analysis::getInstance()
    {
        if (instance_ == nullptr)
        {
            instance_ = new Analysis();
        }
        return static_cast<Analysis *>(instance_);
    }

    void Analysis::selectHadrons()
    {
        selectHadrons_ = true;
    }

    bool Analysis::checkHadron()
    {
        if (int ident = phdst::IPHPIC("IDEN", 0))
        {
            // std::cout << std::hex << phdst::IPILOT(ident+6) << " " << phdst::IPILOT(ident+6) << std::endl;
            // std::cout << std::hex << phdst::IPILOT(ident+6) << " " << phdst::IPILOT(ident+6) << std::endl;
            return (phdst::IPILOT(ident+6) & 0x03) == 0x03;
        } 
        else 
        {
            // poor man hadronic selection
            return NCTR4 >= 5 && ECHAR/ECMAS >= 0.12;
        }
    }

    void Analysis::setOption(const std::string &option, int value)
    {
        if (optionMap_.find(option) != optionMap_.end())
        {
            options_[option] = value;
        }
        else
        {
            std::cerr << "Unknown option: " << option << std::endl;
        }
    }

    void Analysis::setCut(const std::string &cut, float value)
    {
        if (cutMap_.find(cut) != cutMap_.end())
        {
            cuts_[cut] = value;
        }
        else
        {
            std::cerr << "Unknown cut: " << cut << std::endl;
        }
    }

    void Analysis::setIntCut(const std::string &intCut, int value)
    {
        if (intCutMap_.find(intCut) != intCutMap_.end())
        {
            intCuts_[intCut] = value;
        }
        else
        {
            std::cerr << "Unknown intCut: " << intCut << std::endl;
        }
    }

    void Analysis::printSkelanaSettings(std::ostream &os)
    {
        os << "\nSkelana flags:" << std::endl;
        for (const auto &it : optionMap_)
        {
            os << "   " << it.first << " = " << *it.second << std::endl;
        }

        if (IFLSTR > 0)
        {
            os << std::fixed << std::setprecision(3);
            os << "\nTrack selection cuts for this run:" << std::endl;
            os << "   TRKMOM = " << TRKMOM(IFLCUT) << std::endl;
            os << "   TRKMAX = " << TRKMAX(IFLCUT) << std::endl;
            os << "   TRKERR = " << TRKERR(IFLCUT) << std::endl;
            os << "   TRKRPH = " << TRKRPH(IFLCUT) << std::endl;
            os << "   TRKZET = " << TRKZET(IFLCUT) << std::endl;
            os << "   TRKLEN = " << TRKLEN(IFLCUT) << std::endl;
            os << "   TRCCOS = " << TRCCOS(IFLCUT) << std::endl;
            os << "   VDONLY = " << VDONLY(IFLCUT) << std::endl;
            os << "   IDVDWZ = " << IDVDWZ(IFLCUT) << std::endl;
            os << "   IHADRJ = " << IHADRJ(IFLCUT) << std::endl;
            os << "   ISTOEL = " << ISTOEL(IFLCUT) << std::endl;
            os << "   EHPC   = " << EHPC(IFLCUT) << std::endl;
            os << "   EFEMC  = " << EFEMC(IFLCUT) << std::endl;
            os << "   EHAC   = " << EHAC(IFLCUT) << std::endl;
            os << "   ESTIC  = " << ESTIC(IFLCUT) << std::endl;
            os << "   TRNCOS = " << TRNCOS(IFLCUT) << std::endl;
            os << "   RECCAL = " << RECCAL(IFLCUT) << std::endl;
            os << std::defaultfloat;
        }
    }

    void Analysis::user00()
    {
        // std::cout << "Skelana::Analysis::user00: Initialising"  << std::endl;

        // Supress floating point errors
        phdst::PHSET("FPE", 0);
        // Skelana initialization
        PSINI();

        for (const auto &opt : options_)
        {   
            if (auto it = optionMap_.find(opt.first); it != optionMap_.end())
            {
                *it->second = opt.second;
            }
        }
        for (const auto &opt : cuts_)
          {
            if (auto it = cutMap_.find(opt.first); it != cutMap_.end())
              {
                *it->second = opt.second;
              }
          }
        for (const auto &opt : intCuts_)
        {
            if (auto it = intCutMap_.find(opt.first); it != intCutMap_.end())
            {
                *it->second = opt.second;
            }
        }

        printSkelanaSettings(std::cout);

        // Read the energy correction  2000
        if (IFLENR > 0)
        {
            PSENRG();
        }

        if (IFLBHP > 0)
        {
            PSBHPC();
        }

        // RunQuality slection
        if (IFLRNQ > 0)
        {
            //TODO: Should be derived from config
            IRQDMN(1) = 1;
            IRQDMN(2) = 1;
            IRQDMN(5) = 7;
            IRQDMN(6) = 7;
            // Read runquality file
            PSRUNQ(0);
        }


    }

    int Analysis::user01()
    {
        // std::cout << "Skelana::Analysis::user01: Processing pilot " << phdst::NEVENT+1 << std::endl;
        if (int need = super::user01(); need < 1)
        {
            return need;
        }

        if (!checkRecordType())
        {
            return 0;
        }

        // RunSelection
        if (IFLRNQ > 0 && phdst::IIIRUN > 0 && PSRUNS() == 0)
        {
            return 0;
        }

        // Bad HPC in 97
        if (IFLBHP && !checkBadHPC())
        {
            return 0;
        }

        // Hadron Selection
        if(selectHadrons_ && !checkHadron())
        {
            return 0;
        }
        return 1;

    }

    bool Analysis::checkRecordType()
    {
        static const std::set<std::string> skipRecType {"RAW", "TAN", "SOR", "COR", "EOR", "SOS", "EOS","BOF","CPT","0072"};

        const std::string recType = phdst::PHRTY();
        if (skipRecType.find(recType) != skipRecType.end())
        {
            return false;
        }
        if (recType != "DST")
        {
            std::cerr << "Unknown record type: " << recType << std::endl;
            return false;
        }
        return true;
    }

    bool Analysis::checkBadHPC()
    {
        if (phdst::IIIRUN >= IRNBAD(1) && phdst::IIIRUN <= IRNBAD(NEVBAD))
        {
            for (int i = 1; i <= NEVBAD; ++i)
            {
                if (phdst::IIIRUN == IRNBAD(i) && phdst::IIIEVT == IEVBAD(i))
                {
                    return false;
                }
            }
        }
        return true;
    }

    void Analysis::user02()
    {
        // std::cout << "Skelana::Analysis::user02: Processing event " << phdst::NEVENT+1 << std::endl;
        PSBEG();
    }

    void Analysis::user99()
    {
        // std::cout << "Skelana::Analysis::user99: Terminating"  << std::endl;
    }
};
