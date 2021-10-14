/**
 * Copyright 1997-2021 <a href="mailto:zhanghn@zju.edu.cn">Zhang Hongning</a>.
 * 
 * Modified at 2020-02-04
 * Licensed under the Apache License, Version 2.0 (thie "License");
 * You may not use this file except in compliance with the license.
 * You may obtain a copy of the License at
 * 
 *       http://wwww.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language govering permissions and
 * limitations under the License.
 */

package biojar.function.graphics;

import java.awt.Font;
import java.awt.FontFormatException;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.metadata.IIOInvalidTreeException;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageOutputStream;
import org.w3c.dom.Element;
import org.sourceforge.jlibeps.epsgraphics.EpsGraphics2D;
/**
 * 通用图像处理工具，为具体业务提供图像IO等一般性功能。
 * @version 1.0
 * @since 14 2021-02-04
 * @author <a href="mailto:zhanghn@zju.edu.cn">Zhang Hongning</a>
 */
public class ImageHelp {
	/**
	 * 返回特定字体下字符串宽度
	 * @param font 字体
	 * @param string 字符串
	 * @return 宽度
	 */
	public static int getStringFontWidth(Font font, String string) {
		BufferedImage bi = new BufferedImage(1, 1, BufferedImage.TYPE_3BYTE_BGR);
		return bi.createGraphics().getFontMetrics(font).stringWidth(string);
	}
	/**
	 * 将RGB模式的image保存为jpg格式图片
	 * @param image BufferedImage对象
	 * @param filename 输出文件名，建议以jpg/jpeg结尾以便查看器查看
	 * @param dpi 图片水平/垂直分辨率，单位为像素/英寸
	 * @throws IOException 输入输出异常
	 */
	public static void saveAsJPEG(BufferedImage image, String filename, int dpi) throws IOException {
		saveAsJPEG(image, new FileOutputStream(filename), dpi);
	}
	/**
	 * 将RGB模式的image保存为jpg格式图片
	 * @param image 图片对象
	 * @param fos 文件输出流
	 * @param dpi 图像DPI，仅限于jpg
	 * @throws FileNotFoundException 文件缺失异常
	 * @throws IOException 输入输出异常
	 */
	public static void saveAsJPEG(BufferedImage image, FileOutputStream fos, int dpi) throws FileNotFoundException, IOException {
		for (Iterator <ImageWriter> iw = ImageIO.getImageWritersBySuffix("jpg"); iw.hasNext();) {
			ImageWriter writer = iw.next();
			ImageWriteParam writeParam = writer.getDefaultWriteParam();
			writeParam.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
			//调整图片质量
			writeParam.setCompressionQuality(1f);
			IIOMetadata data = writer.getDefaultImageMetadata(new ImageTypeSpecifier(image), writeParam);
			Element tree = (Element) data.getAsTree("javax_imageio_jpeg_image_1.0");
			Element jfif = (Element) tree.getElementsByTagName("app0JFIF").item(0);
			jfif.setAttribute("Xdensity", "" + dpi);
			jfif.setAttribute("Ydensity", "" + dpi);
			jfif.setAttribute("resUnits", "1");// density is dots per inch，如果没有设置会无效
			data.setFromTree("javax_imageio_jpeg_image_1.0", tree);//将tree的内容保存回data，两者无映射关系需此操作，see: http://www.voidcn.com/article/p-zdkeyptk-bts.html
			//输出图片
			ImageOutputStream ios = ImageIO.createImageOutputStream(fos);
			writer.setOutput(ios);
			writer.write(null, new IIOImage(image, null, data), writeParam);
			ios.close();
		}
	}
	/**
	 * 将EpsGraphics2D绘图结果输出为EPS文档，为矢量图格式
	 * @param eps 待输出对象
	 * @param outputfile 输出文件名
	 * @throws IOException 输入输出异常
	 */
	public static void saveAsEPS(EpsGraphics2D eps, String outputfile) throws IOException {
		FileWriter fos = new FileWriter(outputfile);
		fos.write(eps.toString());
		fos.close();
	}
	/**
	 * 将BuffferedImage对象输出为EPS文档，但不改变image本身的位图特征。
	 * @param image 图像对象
	 * @param outputfile 输出文件名
	 * @throws IOException 输入输出异常
	 */
	public static void saveAsEPS(BufferedImage image, String outputfile) throws IOException {
		EpsGraphics2D epsg2d = new EpsGraphics2D();
		epsg2d.drawImage(image, -1, -1, null);
		FileWriter fos = new FileWriter(outputfile);
		fos.write(epsg2d.toString());
		fos.close();
	
	}
	/**
	 * 将RGB模式BufferImage输出文png，并设定dpi
	 * @param image 待输出BufferImage对象
	 * @param fos 文件输出流
	 * @param dpi 输出dpi
	 * @throws IIOInvalidTreeException IIOMeta对象解析异常
	 * @throws IOException 输入输出异常
	 */
	public static void saveAsPNG(BufferedImage image, FileOutputStream fos, int dpi) throws IIOInvalidTreeException, IOException {
		for (Iterator <ImageWriter> iw = ImageIO.getImageWritersByFormatName("png"); iw.hasNext();) {
			ImageWriter writer = iw.next();
			ImageWriteParam writeParam = writer.getDefaultWriteParam();
			ImageTypeSpecifier typeSpecifier = ImageTypeSpecifier.createFromBufferedImageType(BufferedImage.TYPE_INT_RGB);
			IIOMetadata metadata = writer.getDefaultImageMetadata(typeSpecifier, writeParam);
			if (metadata.isReadOnly()||!metadata.isStandardMetadataFormatSupported()) continue;
			
			double inch2cm = 2.54;
			double dotsPerMilli = 1.0 *dpi/ 10 / inch2cm;
			IIOMetadataNode horiz = new IIOMetadataNode("HorizontalPixelSize");  
			horiz.setAttribute("value", Double.toString(dotsPerMilli));
			IIOMetadataNode vert = new IIOMetadataNode("VerticalPixelSize");  
			vert.setAttribute("value", Double.toString(dotsPerMilli));
			IIOMetadataNode dim = new IIOMetadataNode("Dimension");  
			dim.appendChild(horiz);
			dim.appendChild(vert);
			IIOMetadataNode root = new IIOMetadataNode("javax_imageio_1.0");  
			root.appendChild(dim);
			metadata.mergeTree("javax_imageio_1.0", root);
			
			ImageOutputStream ios = ImageIO.createImageOutputStream(fos);
			writer.setOutput(ios);
			writer.write(metadata, new IIOImage(image, null, metadata), writeParam);
			ios.close();
		}
	}
	/**
	 * 获取系统所含所有字体的family名称
	 * @return 代表所有字体family名称的字符串组
	 */
	public static String[] getSystemFont() {
		GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
		return ge.getAvailableFontFamilyNames();
	}
	/**
	 * 获取系统字体路径下所有的ttf字体文件及其对应的字体family组成的哈希表，仅支持Windows和Linux
	 * @return 以字体family名称为Key值，所有相应ttf文件路径名为Value值的哈希表
	 * @throws FontFormatException 字体格式异常
	 * @throws IOException 输入输出异常
	 * @throws Exception 其他异常
	 */
	@SuppressWarnings("unused")
	public static HashMap <String, ArrayList<String>> getSystemFontMap() throws FontFormatException, IOException, Exception {
		String fontdir = null;
		String os = System.getProperty("os.name").toLowerCase();
		if (os.contains("windows")) {
			fontdir = "C:/Windows/Fonts";
		} else if (os.contains("linux")) {
			fontdir = "/usr/share/fonts/";
		} else {
			throw new Exception ("Unsupport Operation System");
		}
		HashMap <String, ArrayList<String>> map = new HashMap<>();
		File[] fontFiles = new File(fontdir).listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().toLowerCase().endsWith("ttf");
			}
		});
		for (File fontfile: fontFiles) {
			String filename = fontfile.getAbsolutePath();
			Font font = Font.createFont(Font.TRUETYPE_FONT, fontfile);
			String fontname = font.getFamily();
			map.putIfAbsent(fontname, new <String> ArrayList<String>());
			map.get(fontname).add(filename);
		
		}
		return map;
	}
	/**
	 * 返回特定字体family下的特定类型字体文件路径名，如Times New Roman, normal输入在Windows下返回
	 * C:/Windows/Fonts/times.ttf，在Linux下返回/usr/share/fonts/times.ttf
	 * @param fontFamilyName 字体family名称
	 * @param fontType 字体类型，分为正常（n、norm、normal），加粗（b、bd、bond），意大利斜体（i、it、Italy）
	 * 和意大利斜体加粗（bi、bondItaly）
	 * @return 获得的字体ttf文件路径名
	 * @throws IOException 输入输出异常
	 * @throws Exception 其他异常
	 */
	public static String getFontFileName(String fontFamilyName, String fontType) throws IOException, Exception {
		HashMap <String, ArrayList<String>> fontMap = getSystemFontMap();
		String normal = null, bond=null, italy=null, bonditaly = null;
		if (!fontMap.containsKey(fontFamilyName)) return null;
		for (String fontfilename: fontMap.get(fontFamilyName)) {//选择加粗字体
			if (fontfilename.toLowerCase().endsWith("b.ttf")||fontfilename.toLowerCase().endsWith("bd.ttf")) 
				bond = fontfilename;
			else if (fontfilename.toLowerCase().endsWith("bi.ttf")) bonditaly = fontfilename;
			else if (fontfilename.toLowerCase().endsWith("i.ttf")) italy = fontfilename;
			else normal = fontfilename;
		}
		switch (fontType) {
			case "b":;
			case "bd":;
			case "bond": return bond;
			case "i":;
			case "it":;
			case "italy": return italy;
			case "bi":;
			case "bonditaly": return bonditaly;
			case "n":;
			case "norm":;
			case "normal":return normal;
			default:return null;
		}
	}
	/*
	 * 将RGB模式的image保存为jpg格式图片
	 * @param image BufferedImage对象
	 * @param outputfile 输出文件对象
	 * @param dpi 图片水平/垂直分辨率，单位为像素/英寸
	 * @throws IOException 输入输出异常
	 * @deprecated 从java 1.7开始，com.sun.image.codec.jpeg被放弃 
	@Deprecated
	public void saveAsJPEG(BufferedImage image, File outputfile, int dpi) throws IOException {
		//import com.sun.image.codec.jpeg.JPEGCodec;
		//import com.sun.image.codec.jpeg.JPEGEncodeParam;
		//import com.sun.image.codec.jpeg.JPEGImageEncoder;
		FileOutputStream fos = new FileOutputStream(outputfile);
		JPEGImageEncoder jpegEncoder = JPEGCodec.createJPEGEncoder(fos);
		JPEGEncodeParam jpegEncodeParam = jpegEncoder.getDefaultJPEGEncodeParam(image);
		jpegEncodeParam.setDensityUnit(JPEGEncodeParam.DENSITY_UNIT_DOTS_INCH);
		jpegEncodeParam.setXDensity(dpi);
		jpegEncodeParam.setYDensity(dpi);
		jpegEncoder.encode(image, jpegEncodeParam);
		fos.close();
	} */
}