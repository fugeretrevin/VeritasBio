import { initializeApp } from "firebase/app";
import { getAuth, GoogleAuthProvider } from "firebase/auth";

// PanoptesBio Web App Firebase configuration
export const firebaseConfig = {
  apiKey: import.meta.env.VITE_FIREBASE_API_KEY || "",
  authDomain: import.meta.env.VITE_FIREBASE_AUTH_DOMAIN || "panoptesbio.firebaseapp.com",
  projectId: import.meta.env.VITE_FIREBASE_PROJECT_ID || "panoptesbio",
  storageBucket: import.meta.env.VITE_FIREBASE_STORAGE_BUCKET || "panoptesbio.firebasestorage.app",
  messagingSenderId: import.meta.env.VITE_FIREBASE_MESSAGING_SENDER_ID || "561444549088",
  appId: import.meta.env.VITE_FIREBASE_APP_ID || "1:561444549088:web:e1be29926bdaaf9766fb85",
  measurementId: import.meta.env.VITE_FIREBASE_MEASUREMENT_ID || "G-3CQM0RGVH3"
};

// Initialize Firebase
const app = initializeApp(firebaseConfig);
export const auth = getAuth(app);
export const googleProvider = new GoogleAuthProvider();

export default app;
