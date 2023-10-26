<!-- ============================================================= -->
<!--                     MATHML 3.0 MODULES                        -->
<!-- ============================================================= -->
<!--                    MATHML 3.0 QUALIFIED NAMES                 -->
<!ENTITY % mathml-qname.mod
                        PUBLIC
"-//W3C//ENTITIES MathML 3.0 Qualified Names 1.0//EN"
"mathml3-qname1.mod"                                                 >


<!--                    MATHML 3.0 DTD                             -->
<!ENTITY % mathml.dtd   PUBLIC
"-//W3C//DTD MathML 3.0//EN" 
"mathml3.dtd"                                                        >


<!--                    MATHML 3.0 EXTRA ENTITIES                  -->
<!ENTITY % ent-mmlextra
                        PUBLIC
"-//W3C//ENTITIES Extra for MathML 3.0//EN" 
"mathml/mmlextra.ent"                                                >


<!--                    MATHML 3.0 ALIASES                         -->
<!ENTITY % ent-mmlalias
                        PUBLIC
"-//W3C//ENTITIES Aliases for MathML 3.0//EN" 
"mathml/mmlalias.ent"                                                >

<!--                    NAMED CHARACTERS                           -->
<!ENTITY % isobox PUBLIC "-//W3C//ENTITIES Box and Line Drawing//EN" "iso8879/isobox.ent">
<!ENTITY % isocyr1 PUBLIC "-//W3C//ENTITIES Russian Cyrillic//EN" "iso8879/isocyr1.ent">
<!ENTITY % isocyr2 PUBLIC "-//W3C//ENTITIES Non-Russian Cyrillic//EN" "iso8879/isocyr2.ent">
<!ENTITY % isodia PUBLIC "-//W3C//ENTITIES Diacritical Marks//EN" "iso8879/isodia.ent">
<!ENTITY % isolat1 PUBLIC "-//W3C//ENTITIES Added Latin 1//EN" "iso8879/isolat1.ent">
<!ENTITY % isolat2 PUBLIC "-//W3C//ENTITIES Added Latin 2//EN" "iso8879/isolat2.ent">
<!ENTITY % isonum PUBLIC "-//W3C//ENTITIES Numeric and Special Graphic//EN" "iso8879/isonum.ent">
<!ENTITY % isopub PUBLIC "-//W3C//ENTITIES Publishing//EN" "iso8879/isopub.ent">
<!ENTITY % isoamsa PUBLIC "-//W3C//ENTITIES Added Math Symbols: Arrow Relations//EN" "iso9573-13/isoamsa.ent">
<!ENTITY % isoamsb PUBLIC "-//W3C//ENTITIES Added Math Symbols: Binary Operators//EN" "iso9573-13/isoamsb.ent">
<!ENTITY % isoamsc PUBLIC "-//W3C//ENTITIES Added Math Symbols: Delimiters//EN" "iso9573-13/isoamsc.ent">
<!ENTITY % isoamsn PUBLIC "-//W3C//ENTITIES Added Math Symbols: Negated Relations//EN" "iso9573-13/isoamsn.ent">
<!ENTITY % isoamso PUBLIC "-//W3C//ENTITIES Added Math Symbols: Ordinary//EN" "iso9573-13/isoamso.ent">
<!ENTITY % isoamsr PUBLIC "-//W3C//ENTITIES Added Math Symbols: Relations//EN" "iso9573-13/isoamsr.ent">
<!ENTITY % isogrk3 PUBLIC "-//W3C//ENTITIES Greek Symbols//EN" "iso9573-13/isogrk3.ent">
<!ENTITY % isomfrk PUBLIC "-//W3C//ENTITIES Math Alphabets: Fraktur//EN" "iso9573-13/isomfrk.ent">
<!ENTITY % isomopf PUBLIC "-//W3C//ENTITIES Math Alphabets: Open Face//EN" "iso9573-13/isomopf.ent">
<!ENTITY % isomscr PUBLIC "-//W3C//ENTITIES Math Alphabets: Script//EN" "iso9573-13/isomscr.ent">
<!ENTITY % isotech PUBLIC "-//W3C//ENTITIES General Technical//EN" "iso9573-13/isotech.ent">


<!-- ============================================================= -->
<!--                    SET UP "mml" AS THE MATH PREFIX            -->
<!-- ============================================================= -->


<!--                    MAKE MATH PREFIX PARAMETER ENTITY HAPPEN   -->
<!ENTITY % MATHML.prefixed 
                        "INCLUDE"                                    >


<!--                    MAKE PREFIX EQUAL "mml"                    -->
<!ENTITY % MATHML.prefix "mml"                                       >


<!--                    SET UP "pfx"                               -->
<![%MATHML.prefixed;[
<!ENTITY % MATHML.pfx   "%MATHML.prefix;:"                           >
]]>


<!--                    USE "pfx" TO SET THE MATH ELEMENT NAME     -->
<!ENTITY % math.qname   "%MATHML.pfx;math"                           >



<!-- ============================================================= -->
<!--                    CALL THE MATHML ENTITIES                   -->
<!-- ============================================================= -->


<!--                    MATHML CHARACTER ENTITIES                  -->
<!--                    Set the "INCLUDE" or "IGNORE" marked section
                        switch for the MATHML 3.0 DTD. This switch
                        determines whether the math processing in
                        the -%mathml.dtd; module or the
                        -%xmlspecchars.ent; module in this DTD Suite
                        will invoke the sets of special characters
                        that are defined AND INVOKED in the
                        -%xmlspecchars.ent; module.
                        A value of "IGNORE" turns off the second
                        invocation from the MathML DTD.            -->
<!ENTITY % mathml-charent.module
                        "IGNORE"                                     >


<!--                    MATHML-SPECIFIC CHARACTER ENTITIES         -->
<!--                    Because the MATHML invocation is canceled
                        using the parameter entity just defined,
                        the two external entities below must be
                        invoked from here, as neither XMLspecchars
                        nor the MathML DTD will invoke them.       -->

<!--                    New characters defined by MathML           -->
%ent-mmlextra;


<!--                    MathML aliases for characters defined above-->
%ent-mmlalias;

<!--                    NAMED CHARACTERS                           -->
%isobox;
%isocyr1;
%isocyr2;
%isodia;
%isolat1;
%isolat2;
%isonum;
%isopub;
%isoamsa;
%isoamsb;
%isoamsc;
%isoamsn;
%isoamso;
%isoamsr;
%isogrk3;
%isomfrk;
%isomopf;
%isomscr;
%isotech;

<!-- ============================================================= -->
<!--                    RESTRICT CONTENT OF ANNOTATION             -->
<!-- ============================================================= -->

<!--                    MATHML ANNOTATION MODEL                    -->
<!--                    The MathML 3.0 DTD establishes the content of 
                        the <mml:annotation-xml> element as any MathML
                        expression. This is unnecessarily broad, 
                        causing inconvenience during conversion. 
                        The model has been restricted to
                        one or more paragraphs.                    -->
<!ENTITY % annotation-xml.model
                        "#PCDATA"                                       >


<!-- ============================================================= -->
<!--                    MathML 3.0 INVOCATION                      -->
<!-- ============================================================= -->

%mathml.dtd;
