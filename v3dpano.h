#ifndef V3DPANO_H
#define V3DPANO_H

class v3Dpano
{
  public:
    //v3Dpano(gpano_t *self);
    v3Dpano(int wIm, int hIm, int bIm, int wRa, int hRa, int bRa, int wPl, int hPl, int bPl, int nPano);
    virtual ~v3Dpano();
    
    bool creatPanoBySingle(gpano_t &self, int i);
    bool fillPanoHoleBySingle(gpano_t &self, int i);
    bool newRaIm2PTS(gpano_t &self, string filename);
    
  private:
    int m_nPano;
    int m_stepIm, m_bandIm, m_wIm, m_hIm;
    int m_stepRa, m_bandRa, m_wRa, m_hRa;
    int m_stepPl, m_bandPl, m_wPl, m_hPl;

   uchar *m_imNew;
   ushort *m_raNew;
   
   int *m_imProj[10];
   double *m_curRaRate;
   uchar *m_curImSrc;
   
   bool calNormalVector(ushort *ra, int x, int y, vec3_t &P,vec3_t &v, ushort &r);
};

#endif // V3DPANO_H
